/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "viewFactorSky.H"
#include "surfaceFields.H"
#include "constants.H"
#include "greyDiffusiveViewFactorFixedValueFvPatchScalarField.H"
#include "typeInfo.H"
#include "addToRunTimeSelectionTable.H"

#include "wallFvPatch.H"
#include "TableFile.H"

#include "mappedPatchBase.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
    defineTypeNameAndDebug(viewFactorSky, 0);
    addToRadiationRunTimeSelectionTables(viewFactorSky);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::radiationModels::viewFactorSky::initialise()
{
    const polyBoundaryMesh& coarsePatches = coarseMesh_.boundaryMesh();
    const volScalarField::Boundary& qrp = qr_.boundaryField();

    label count = 0;
    forAll(qrp, patchi)
    {
        // const polyPatch& pp = mesh_.boundaryMesh()[patchi];
        const fvPatchScalarField& qrPatchi = qrp[patchi];

        if ((isA<fixedValueFvPatchScalarField>(qrPatchi)))
        {
            selectedPatches_[count] = qrPatchi.patch().index();
            nLocalCoarseFaces_ += coarsePatches[patchi].size();
            count++;
        }
    }

    selectedPatches_.resize(count--);

    if (debug)
    {
        Pout<< "radiationModels::viewFactorSky::initialise() "
            << "Selected patches:" << selectedPatches_ << endl;
        Pout<< "radiationModels::viewFactorSky::initialise() "
            << "Number of coarse faces:" << nLocalCoarseFaces_ << endl;
    }

    totalNCoarseFaces_ = nLocalCoarseFaces_;
    reduce(totalNCoarseFaces_, sumOp<label>());

    if (debug && Pstream::master())
    {
        InfoInFunction
            << "Total number of clusters : " << totalNCoarseFaces_ << endl;
    }

    labelListIOList subMap
    (
        IOobject
        (
            "subMap",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    labelListIOList constructMap
    (
        IOobject
        (
            "constructMap",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    IOList<label> consMapDim
    (
        IOobject
        (
            "constructMapDim",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    map_.reset
    (
        new mapDistribute
        (
            consMapDim[0],
            move(subMap),
            move(constructMap)
        )
    );

    // Read GMRES parameters (with defaults for backward compatibility)
    gmresRestart_ = coeffs_.lookupOrDefault<label> ("gmresRestart",  20);
    gmresMaxIter_ = coeffs_.lookupOrDefault<label> ("gmresMaxIter",  200);
    gmresTol_     = coeffs_.lookupOrDefault<scalar>("gmresTolerance", 1e-6);

    const bool smoothing = readBool(coeffs_.lookup("smoothing"));

    // Read F (view factor values) and globalFaceFaces (global column indices)
    // written by viewFactorsGen — same files used by directAndDiffuse.
    {
        scalarListIOList F
        (
            IOobject
            (
                "F",
                mesh_.facesInstance(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        labelListIOList globalFaceFaces
        (
            IOobject
            (
                "globalFaceFaces",
                mesh_.facesInstance(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        if (F.size() != nLocalCoarseFaces_)
        {
            FatalErrorInFunction
                << "F size=" << F.size()
                << " != nLocalCoarseFaces=" << nLocalCoarseFaces_
                << ". Re-run viewFactorsGen."
                << exit(FatalError);
        }

        globalIndex globalNumbering(nLocalCoarseFaces_);
        globalOffset_ = globalNumbering.toGlobal(Pstream::myProcNo(), 0);

        FsparseCols_.setSize(nLocalCoarseFaces_);
        FsparseVals_.setSize(nLocalCoarseFaces_);

        label nnzTotal = 0;
        forAll(F, localI)
        {
            FsparseCols_[localI] = globalFaceFaces[localI];
            FsparseVals_[localI] = F[localI];

            if (smoothing)
            {
                scalarList& vals = FsparseVals_[localI];
                scalar sumF = 0.0;
                forAll(vals, k) sumF += vals[k];
                if (sumF > VSMALL)
                {
                    const scalar fac = 1.0 - (sumF - 1.0)/(sumF + 0.001);
                    forAll(vals, k) vals[k] *= fac;
                }
            }
            nnzTotal += FsparseCols_[localI].size();
        }
        Pout << "Read F+globalFaceFaces: " << nLocalCoarseFaces_
             << " rows, nnz=" << nnzTotal
             << ", globalOffset=" << globalOffset_ << endl;
    }

    qPrev_.setSize(totalNCoarseFaces_, 0.0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::viewFactorSky::viewFactorSky(const volScalarField& T)
:
    radiationModel(typeName, T),
    finalAgglom_
    (
        IOobject
        (
            "finalAgglom",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    map_(),
    coarseMesh_
    (
        IOobject
        (
            mesh_.name(),
            mesh_.polyMesh::instance(),
            mesh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        finalAgglom_
    ),
    qr_
    (
        IOobject
        (
            "qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    selectedPatches_(mesh_.boundary().size(), -1),
    totalNCoarseFaces_(0),
    nLocalCoarseFaces_(0),
    globalOffset_(0),
    gmresRestart_(20),
    gmresMaxIter_(200),
    gmresTol_(1e-6)
{
    initialise();
}


Foam::radiationModels::viewFactorSky::viewFactorSky
(
    const dictionary& dict,
    const volScalarField& T
)
:
    radiationModel(typeName, dict, T),
    finalAgglom_
    (
        IOobject
        (
            "finalAgglom",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    map_(),
    coarseMesh_
    (
        IOobject
        (
            mesh_.name(),
            mesh_.polyMesh::instance(),
            mesh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        finalAgglom_
    ),
    qr_
    (
        IOobject
        (
            "qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    selectedPatches_(mesh_.boundary().size(), -1),
    totalNCoarseFaces_(0),
    nLocalCoarseFaces_(0),
    globalOffset_(0),
    gmresRestart_(20),
    gmresMaxIter_(200),
    gmresTol_(1e-6)
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::viewFactorSky::~viewFactorSky()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiationModels::viewFactorSky::read()
{
    if (radiationModel::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}



Foam::scalarField Foam::radiationModels::viewFactorSky::matvecC
(
    const scalarField& x,
    const scalarField& E
) const
{
    // Each proc computes w for its local rows [globalOffset_, globalOffset_+nLocal).
    // AllGather (sumOp) assembles the full w on all procs.
    scalarField w(totalNCoarseFaces_, 0.0);

    for (label localI = 0; localI < nLocalCoarseFaces_; localI++)
    {
        const label globalI = globalOffset_ + localI;
        const labelList& cols = FsparseCols_[localI];
        const scalarList& vals = FsparseVals_[localI];

        scalar Fii = 0.0;
        scalar offdiag = 0.0;
        forAll(cols, k)
        {
            const label j = cols[k];
            const scalar Fij = vals[k];
            if (j == globalI)
            {
                Fii = Fij;
            }
            else
            {
                offdiag += (1.0 - 1.0/E[j])*Fij*x[j];
            }
        }
        const scalar invEi = 1.0/E[globalI];
        w[globalI] = (invEi - (invEi - 1.0)*Fii)*x[globalI] + offdiag;
    }

    // AllGather: each proc filled only its chunk (rest is 0), sum assembles full w
    Pstream::listCombineGather(w, plusEqOp<scalar>());  // sumOp non modifica il primo argomento: la raccolta non sommava nulla
    Pstream::listCombineScatter(w);

    return w;
}


Foam::scalar Foam::radiationModels::viewFactorSky::solveGMRES
(
    const scalarField& b,
    scalarField& x,
    const scalarField& E
) const
{
    const label N = b.size();
    const label m = min(gmresRestart_, N);

    // Absolute tolerance derived from initial RHS norm
    const scalar bnorm = Foam::sqrt(sum(b*b));
    const scalar tolAbs = gmresTol_ * max(bnorm, VSMALL);

    scalarField r = b - matvecC(x, E);
    scalar beta = Foam::sqrt(sum(r*r));

    if (beta <= tolAbs)
        return beta;

    // Krylov basis vectors V[0..m]
    List<scalarField> V(m + 1, scalarField(N, 0.0));
    // Upper Hessenberg H stored row-major: H[row][col], size (m+1) x m
    List<scalarField> H(m + 1, scalarField(m, 0.0));
    // Givens rotation coefficients and RHS vector g
    scalarField cs(m, 0.0), sn(m, 0.0), g(m + 1, 0.0);

    label totalIter = 0;

    for (label outer = 0; outer < gmresMaxIter_; outer++)
    {
        // (Re-)initialize Krylov space from current residual
        forAll(V, i) V[i] = scalar(0.0);
        forAll(H, i) H[i] = scalar(0.0);
        cs = scalar(0.0);
        sn = scalar(0.0);
        g  = scalar(0.0);

        V[0] = r / beta;
        g[0] = beta;

        label j = 0;
        for (; j < m; j++)
        {
            scalarField w = matvecC(V[j], E);

            // Modified Gram-Schmidt orthogonalization
            for (label i = 0; i <= j; i++)
            {
                H[i][j] = sum(w*V[i]);
                w -= H[i][j]*V[i];
            }
            H[j+1][j] = Foam::sqrt(sum(w*w));

            if (H[j+1][j] > VSMALL)
                V[j+1] = w/H[j+1][j];

            // Apply previous Givens rotations to new column j of H
            for (label i = 0; i < j; i++)
            {
                const scalar tmp =  cs[i]*H[i][j] + sn[i]*H[i+1][j];
                H[i+1][j]       = -sn[i]*H[i][j] + cs[i]*H[i+1][j];
                H[i][j]         = tmp;
            }

            // Compute and apply new Givens rotation to zero H[j+1][j]
            const scalar denom =
                Foam::sqrt(H[j][j]*H[j][j] + H[j+1][j]*H[j+1][j]);
            if (denom > VSMALL)
            {
                cs[j] = H[j][j]/denom;
                sn[j] = H[j+1][j]/denom;
            }
            else
            {
                cs[j] = 1.0; sn[j] = 0.0;
            }
            H[j][j]   =  cs[j]*H[j][j] + sn[j]*H[j+1][j];
            H[j+1][j] = 0.0;

            g[j+1] = -sn[j]*g[j];
            g[j]   =  cs[j]*g[j];

            totalIter++;
            beta = mag(g[j+1]);

            if (beta <= tolAbs)
            {
                j++;
                break;
            }
        }

        // Back-substitute upper-triangular system H[0..j-1][0..j-1] * y = g[0..j-1]
        scalarField y(j, 0.0);
        for (label i = j - 1; i >= 0; i--)
        {
            y[i] = g[i];
            for (label k = i + 1; k < j; k++)
                y[i] -= H[i][k]*y[k];
            if (mag(H[i][i]) > VSMALL)
                y[i] /= H[i][i];
        }

        // Update solution: x = x + V[0..j-1] * y
        for (label i = 0; i < j; i++)
            x += y[i]*V[i];

        if (beta <= tolAbs)
            break;

        // Restart: recompute residual
        r = b - matvecC(x, E);
        beta = Foam::sqrt(sum(r*r));

        if (beta <= tolAbs)
            break;
    }

    if (Pstream::master())
    {
        Info<< "    GMRES: " << totalIter << " iterations, |r| = " << beta << endl;
    }
    return beta;
}


void Foam::radiationModels::viewFactorSky::calculate()
{
    // Store previous iteration
    qr_.storePrevIter();

    scalarField compactCoarseT4(map_->constructSize(), 0.0);
    scalarField compactCoarseE(map_->constructSize(), 0.0);
    scalarField compactCoarseHo(map_->constructSize(), 0.0);

    globalIndex globalNumbering(nLocalCoarseFaces_);

    // Fill local averaged(T), emissivity(E) and external heatFlux(Ho)
    DynamicList<scalar> localCoarseT4ave(nLocalCoarseFaces_);
    DynamicList<scalar> localCoarseEave(nLocalCoarseFaces_);
    DynamicList<scalar> localCoarseHoave(nLocalCoarseFaces_);

    volScalarField::Boundary& qrBf = qr_.boundaryFieldRef();
    
    //////////////////////////////////////////////////////////////////////////
    //obtain Tambient to calculate Tsky - can find a better way to import Tambient?
    Time& time = const_cast<Time&>(mesh_.time());
    //label timestep = ceil( (time.value()/3600)-1E-6 ); timestep = timestep%24;
    
    dictionary TambientIO;
    TambientIO.add(
        "file", 
        fileName
        (
            "$FOAM_CASE/0/air/Tambient"
        )
    );
    Function1s::TableFile<scalar> Tambient
    (
        "Tambient",
        TambientIO
    );
    //////////////////////////////////////////////////////////////////////////
    fileName cloudCoverFile
    (
       "$FOAM_CASE/0/air/cloudCover"
    );
    
    scalar cc = 0; //cloud cover
    if(isFile(cloudCoverFile))
    {
        Info << "Reading cloud cover values..." << endl;
        dictionary cloudCoverIO;
        cloudCoverIO.add(
            "file", 
            cloudCoverFile
        );
        Function1s::TableFile<scalar> cloudCover
        (
            "cloudCover",
            cloudCoverIO
        );        
        cc = cloudCover.value(time.value());
    }
    else
    {
        Info << "Constant cloud cover of 0 is being used..." << endl;
    }
    
    //////////////////////////////////////////////////////////////////////////
    //is grass model activated?
    const polyMesh& airMesh =
        mesh_.time().lookupObject<polyMesh>("air");
    IOdictionary grassProperties
    (
        IOobject
        (
            "grassProperties",
            airMesh.time().constant(),
            airMesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );
    if (grassProperties.typeHeaderOk<IOdictionary>(true))
    {
        word grassModel(grassProperties.lookup("grassModel"));
        if (grassModel != "none")
        {
            const dictionary& modelCoeffs = grassProperties.subDict(grassModel + "Coeffs");
            const List<word>& grassPatches_ = modelCoeffs.lookup("grassPatches");
            grassPatches = grassPatches_;
        }
    }
    //////////////////////////////////////////////////////////////////////////

    forAll(selectedPatches_, i)
    {
        label patchID = selectedPatches_[i];

        const scalarField& Tp = T_.boundaryField()[patchID];
        const scalarField& sf = mesh_.magSf().boundaryField()[patchID];
        scalarField Tg_(T_.size(), -1.0);
        if (grassPatches.found(mesh_.boundary()[patchID].name()))
        {
            if(mesh_.name() == "vegetation")
            {
                // Get the coupling information from the mappedPatchBase
                const mappedPatchBase& mpp =
                    refCast<const mappedPatchBase>(mesh_.boundary()[patchID].patch());
                const polyMesh& nbrMesh = mpp.sampleMesh();
                const label samplePatchI = mpp.samplePolyPatch().index();
                const fvPatch& nbrPatch =
                    refCast<const fvMesh>(nbrMesh).boundary()[samplePatchI];
                scalarField TgNbr = nbrPatch.lookupPatchField<volScalarField, scalar>("Tg");
                    mpp.distribute(TgNbr);
                Tg_ = TgNbr;
            }
            else
            {
                Tg_ = mesh_.boundary()[patchID].lookupPatchField<volScalarField, scalar>("Tg");
            }
        }

        fvPatchScalarField& qrPatch = qrBf[patchID];

        greyDiffusiveViewFactorFixedValueFvPatchScalarField& qrp =
            refCast
            <
                greyDiffusiveViewFactorFixedValueFvPatchScalarField
            >(qrPatch);

        const scalarList eb = qrp.emissivity();

        const scalarList& Hoi = qrp.qro();

        const polyPatch& pp = coarseMesh_.boundaryMesh()[patchID];
        const labelList& coarsePatchFace = coarseMesh_.patchFaceMap()[patchID];

        scalarList T4ave(pp.size(), 0.0);
        scalarList Eave(pp.size(), 0.0);
        scalarList Hoiave(pp.size(), 0.0);

        if (pp.size() > 0)
        {
            const labelList& agglom = finalAgglom_[patchID];
            label nAgglom = max(agglom) + 1;

            labelListList coarseToFine(invertOneToMany(nAgglom, agglom));

            forAll(coarseToFine, coarseI)
            {
                const label coarseFaceID = coarsePatchFace[coarseI];
                const labelList& fineFaces = coarseToFine[coarseFaceID];
                UIndirectList<scalar> fineSf
                (
                    sf,
                    fineFaces
                );

                const scalar area = sum(fineSf());

                // Temperature, emissivity and external flux area weighting
                forAll(fineFaces, j)
                {
                    label facei = fineFaces[j];
                    if (!isA<wallFvPatch>(mesh_.boundary()[patchID])) // added to take into account sky temperature
                    {
                        scalar Tambient_ = Tambient.value(time.value());
                        scalar ec = (1-0.84*cc)*(0.527 + 0.161*Foam::exp(8.45*(1-273/Tambient_))) +0.84*cc; //cloud emissivity
                        scalar Tsky = pow(9.365574E-6*(1-cc)*pow(Tambient_,6) + pow(Tambient_,4)*cc*ec ,0.25); // Swinbank model (1963, Cole 1976)
                        T4ave[coarseI] += (pow4(Tsky)*sf[facei])/area;
                    }
                    else
                    {
                        if (grassPatches.found(mesh_.boundary()[patchID].name()))//use Tg if patch is covered with grass
                        {
                            T4ave[coarseI] += (pow4(Tg_[facei])*sf[facei])/area;
                        }
                        else//otherwise use T wall temperature
                        {
                            T4ave[coarseI] += (pow4(Tp[facei])*sf[facei])/area;
                        }
                    }
                    Eave[coarseI] += (eb[facei]*sf[facei])/area;
                    Hoiave[coarseI] += (Hoi[facei]*sf[facei])/area;
                }
            }
        }

        localCoarseT4ave.append(T4ave);
        localCoarseEave.append(Eave);
        localCoarseHoave.append(Hoiave);
    }

    // Fill the local values to distribute
    SubList<scalar>(compactCoarseT4, nLocalCoarseFaces_) = localCoarseT4ave;
    SubList<scalar>(compactCoarseE, nLocalCoarseFaces_) = localCoarseEave;
    SubList<scalar>(compactCoarseHo, nLocalCoarseFaces_) = localCoarseHoave;

    // Distribute data
    map_->distribute(compactCoarseT4);
    map_->distribute(compactCoarseE);
    map_->distribute(compactCoarseHo);

    // Distribute local global ID
    labelList compactGlobalIds(map_->constructSize(), 0.0);

    labelList localGlobalIds(nLocalCoarseFaces_);

    for(label k = 0; k < nLocalCoarseFaces_; k++)
    {
        localGlobalIds[k] = globalNumbering.toGlobal(Pstream::myProcNo(), k);
    }

    SubList<label>
    (
        compactGlobalIds,
        nLocalCoarseFaces_
    ) = localGlobalIds;

    map_->distribute(compactGlobalIds);

    // Create global size vectors
    scalarField T4(totalNCoarseFaces_, 0.0);
    scalarField E(totalNCoarseFaces_, 0.0);
    scalarField qrExt(totalNCoarseFaces_, 0.0);

    // Fill lists from compact to global indexes.
    forAll(compactCoarseT4, i)
    {
        T4[compactGlobalIds[i]] = compactCoarseT4[i];
        E[compactGlobalIds[i]] = compactCoarseE[i];
        qrExt[compactGlobalIds[i]] = compactCoarseHo[i];
    }

    Pstream::listCombineGather(T4, maxEqOp<scalar>());
    Pstream::listCombineGather(E, maxEqOp<scalar>());
    Pstream::listCombineGather(qrExt, maxEqOp<scalar>());

    Pstream::listCombineScatter(T4);
    Pstream::listCombineScatter(E);
    Pstream::listCombineScatter(qrExt);

    // Build RHS: each proc fills its local rows, then AllGather
    // b[i] = sum_j F[i][j]*sigma*T4[j] - sigma*T4[i] - qrExt[i]
    const scalar sigma = physicoChemical::sigma.value();
    scalarField b(totalNCoarseFaces_, 0.0);
    for (label localI = 0; localI < nLocalCoarseFaces_; localI++)
    {
        const label globalI = globalOffset_ + localI;
        b[globalI] = -sigma*T4[globalI] - qrExt[globalI];
        const labelList& cols = FsparseCols_[localI];
        const scalarList& vals = FsparseVals_[localI];
        forAll(cols, k)
        {
            b[globalI] += vals[k]*sigma*T4[cols[k]];
        }
    }
    Pstream::listCombineGather(b, plusEqOp<scalar>());
    Pstream::listCombineScatter(b);

    // Distributed GMRES: all procs participate, warm-start from previous solution
    scalarField x(qPrev_);
    if (Pstream::master())
        Info<< "\nSolving LW view factor equations (GMRES distributed)..." << endl;
    const scalar finalRes = solveGMRES(b, x, E);
    if (Pstream::master())
        Info<< "GMRES converged: residual = " << finalRes << endl;

    // x is full and replicated on all procs — no scatter needed
    scalarField q(x);
    qPrev_ = x;

    label globCoarseId = 0;
    forAll(selectedPatches_, i)
    {
        const label patchID = selectedPatches_[i];
        const polyPatch& pp = mesh_.boundaryMesh()[patchID];
        if (pp.size() > 0)
        {
            scalarField& qrp = qrBf[patchID];
            const scalarField& sf = mesh_.magSf().boundaryField()[patchID];
            const labelList& agglom = finalAgglom_[patchID];
            label nAgglom = max(agglom)+1;

            labelListList coarseToFine(invertOneToMany(nAgglom, agglom));

            const labelList& coarsePatchFace =
                coarseMesh_.patchFaceMap()[patchID];

            // Correzione per faccia fine (come per l'onda corta): il sistema usa per la coarse la media
            // pesata per area di T^4, quindi una faccia fine piu' calda della sua coarse non emetteva
            // l'eccesso (parete a 423 K con qr -3 W/m2). Con irradiazione uniforme sulla coarse:
            // q_fine = q_coarse - eps*sigma*(Tf^4 - <T^4>). La somma pesata per area delle correzioni e' 0.
            const bool fineCorr =
                isA<wallFvPatch>(mesh_.boundary()[patchID])
             && !grassPatches.found(mesh_.boundary()[patchID].name());
            const scalarField& Tp = T_.boundaryField()[patchID];
            const scalarList eb =
                refCast<greyDiffusiveViewFactorFixedValueFvPatchScalarField>
                (
                    qrBf[patchID]
                ).emissivity();

            scalar heatFlux = 0.0;
            forAll(coarseToFine, coarseI)
            {
                label globalCoarse =
                    globalNumbering.toGlobal(Pstream::myProcNo(), globCoarseId);
                const label coarseFaceID = coarsePatchFace[coarseI];
                const labelList& fineFaces = coarseToFine[coarseFaceID];

                scalar T4c = 0;
                if (fineCorr)
                {
                    scalar area = 0;
                    forAll(fineFaces, k)
                    {
                        T4c += pow4(Tp[fineFaces[k]])*sf[fineFaces[k]];
                        area += sf[fineFaces[k]];
                    }
                    T4c /= max(area, vSmall);
                }

                forAll(fineFaces, k)
                {
                    label facei = fineFaces[k];

                    qrp[facei] = q[globalCoarse];
                    if (fineCorr)
                    {
                        qrp[facei] -= eb[facei]*sigma*(pow4(Tp[facei]) - T4c);
                    }
                    heatFlux += qrp[facei]*sf[facei];
                }
                globCoarseId ++;
            }
        }
    }

    if (debug)
    {
        forAll(qrBf, patchID)
        {
            const scalarField& qrp = qrBf[patchID];
            const scalarField& magSf = mesh_.magSf().boundaryField()[patchID];
            const scalar heatFlux = gSum(qrp*magSf);

            InfoInFunction
                << "Total heat transfer rate at patch: "
                << patchID << " "
                << heatFlux << endl;
        }
    }

    // Relax qr if necessary
    qr_.relax();
}


Foam::tmp<Foam::volScalarField> Foam::radiationModels::viewFactorSky::Rp() const
{
    return volScalarField::New
    (
        "Rp",
        mesh_,
        dimensionedScalar
        (
            dimMass/pow3(dimTime)/dimLength/pow4(dimTemperature),
            0
        )
    );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::radiationModels::viewFactorSky::Ru() const
{
    return volScalarField::Internal::New
    (
        "Ru",
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), 0)
    );
}

// ************************************************************************* //
