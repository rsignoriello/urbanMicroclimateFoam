/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "directAndDiffuse.H"
#include "surfaceFields.H"
#include "constants.H"
#include "solarLoadViewFactorFixedValueFvPatchScalarField.H"
#include "wallFvPatch.H"
#include "typeInfo.H"
#include "Time.H"

#include "vectorIOList.H"

#include "TableFile.H"

#include "mappedPatchBase.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace solarLoad
    {
        defineTypeNameAndDebug(directAndDiffuse, 0);
        addToSolarLoadRunTimeSelectionTables(directAndDiffuse);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solarLoad::directAndDiffuse::initialise()
{
    const polyBoundaryMesh& coarsePatches = coarseMesh_.boundaryMesh();
    const volScalarField::Boundary& qsp = qs_.boundaryField();

    label count = 0;
    forAll(qsp, patchI)
    {
        //const polyPatch& pp = mesh_.boundaryMesh()[patchI];
        const fvPatchScalarField& qsPatchI = qsp[patchI];

        if ((isA<fixedValueFvPatchScalarField>(qsPatchI)))
        {
            selectedPatches_[count] = qsPatchI.patch().index();
            nLocalCoarseFaces_ += coarsePatches[patchI].size();

            if ((isA<wallFvPatch>(mesh_.boundary()[patchI])))
            {
                wallPatchOrNot_[count] = 1;
                nLocalWallCoarseFaces_ += coarsePatches[patchI].size();
                nLocalFineFaces_ += qsPatchI.patch().size();
            }

            count++;
        }
    }
    Info<< "Selected patches:" << selectedPatches_ << endl;
    Info<< "Number of coarse faces:" << nLocalCoarseFaces_ << endl;
    Info << "wallPatchOrNot_: " << wallPatchOrNot_ << endl;
    Info << "nLocalWallCoarseFaces_: " << nLocalWallCoarseFaces_ << endl;

    selectedPatches_.resize(count);
    wallPatchOrNot_.resize(count);

    Info<< "Selected patches:" << selectedPatches_ << endl;
    Info<< "Number of coarse faces:" << nLocalCoarseFaces_ << endl;
    Info << "wallPatchOrNot_: " << wallPatchOrNot_ << endl;
    Info << "nLocalWallCoarseFaces_: " << nLocalWallCoarseFaces_ << endl;

    if (debug)
    {
        Pout<< "Selected patches:" << selectedPatches_ << endl;
        Pout<< "Number of coarse faces:" << nLocalCoarseFaces_ << endl;
    }

    totalNCoarseFaces_ = nLocalCoarseFaces_;
    reduce(totalNCoarseFaces_, sumOp<label>());
    totalNFineFaces_ = nLocalFineFaces_;
    reduce(totalNFineFaces_, sumOp<label>());

    if (Pstream::master())
    {
        Info<< "Total number of clusters : " << totalNCoarseFaces_ << endl;
        Info<< "Total number of fine faces : " << totalNFineFaces_ << endl;
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

    scalarListIOList FmyProc
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

    scalarListIOList solarLoadFineFacesmyProc
    (
        IOobject
        (
            "solarLoadFineFaces",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    solarLoadFineFacesSize = solarLoadFineFacesmyProc.size();

    scalarListIOList skyViewCoeffmyProc
    (
        IOobject
        (
            "skyViewCoeff",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    skyViewCoeffSize = skyViewCoeffmyProc.size();

    scalarListIOList sunViewCoeffmyProc
    (
        IOobject
        (
            "sunViewCoeff",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    sunViewCoeffSize = sunViewCoeffmyProc.size();

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

    // Store local view factor data for the parallel Jacobi solver
    FmyProc_ = FmyProc;
    globalFaceFacesLocal_ = globalFaceFaces;

    // Direct local assignment — each proc keeps only its own chunk.
    // No gatherList needed: skyViewCoeff/sunViewCoeff/solarLoadFineFaces are
    // per-face properties; calculate() only accesses local indices.
    solarLoadFineFacesGlobal_.reset(new scalarListList(solarLoadFineFacesSize));
    forAll(solarLoadFineFacesGlobal_(), lo)
        solarLoadFineFacesGlobal_()[lo] = solarLoadFineFacesmyProc[lo];

    skyViewCoeffGlobal_.reset(new scalarListList(skyViewCoeffSize));
    sunViewCoeffGlobal_.reset(new scalarListList(sunViewCoeffSize));
    forAll(skyViewCoeffGlobal_(), lo)
    {
        skyViewCoeffGlobal_()[lo] = skyViewCoeffmyProc[lo];
        sunViewCoeffGlobal_()[lo] = sunViewCoeffmyProc[lo];
    }

    // Apply optional view factor row normalization (smoothing) locally on each
    // proc. Equivalent to global row normalization (scale = 1/rowSum) but avoids
    // allocating the O(N^2) global Fmatrix_ (which would require ~135 GB for
    // N=130k). Each proc normalises only its own rows independently.
    constAlbedo_ = readBool(coeffs_.lookup("constantAlbedo"));
    bool smoothing = readBool(coeffs_.lookup("smoothing"));
    if (smoothing)
    {
        Info<< "Smoothing the view factor rows (local normalization)..." << endl;
        forAll(FmyProc_, localI)
        {
            scalar rowSum = 0.0;
            forAll(FmyProc_[localI], k)
                rowSum += FmyProc_[localI][k];
            if (rowSum > SMALL)
            {
                forAll(FmyProc_[localI], k)
                    FmyProc_[localI][k] /= rowSum;
            }
        }
    }

    // Read Jacobi solver parameters (fall back to defaults if absent)
    maxJacobiIter_ = coeffs_.lookupOrDefault<label>("maxJacobiIter", 200);
    jacobiTol_     = coeffs_.lookupOrDefault<scalar>("jacobiTol", scalar(1.0e-6));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solarLoad::directAndDiffuse::directAndDiffuse(const volScalarField& T)
:
    solarLoadModel(typeName, T),
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
    qs_
    (
        IOobject
        (
            "qs",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    Fmatrix_(),
    CLU_(),
    solarLoadFineFacesGlobal_(),
    skyViewCoeffGlobal_(),
    sunViewCoeffGlobal_(),
    selectedPatches_(mesh_.boundary().size(), -1),
    wallPatchOrNot_(mesh_.boundary().size(), 0),
    totalNCoarseFaces_(0),
    nLocalCoarseFaces_(0),
    nLocalWallCoarseFaces_(0),
    nLocalFineFaces_(0),
    totalNFineFaces_(0),
    constAlbedo_(false),
    timestepsInADay_(24),
    iterCounter_(0),
    pivotIndices_(0),
    FmyProc_(),
    globalFaceFacesLocal_(),
    maxJacobiIter_(200),
    jacobiTol_(1.0e-6)
{
    initialise();
}


Foam::solarLoad::directAndDiffuse::directAndDiffuse
(
    const dictionary& dict,
    const volScalarField& T
)
:
    solarLoadModel(typeName, dict, T),
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
    qs_
    (
        IOobject
        (
            "Qs",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    Fmatrix_(),
    CLU_(),
    solarLoadFineFacesGlobal_(),
    skyViewCoeffGlobal_(),
    sunViewCoeffGlobal_(),
    selectedPatches_(mesh_.boundary().size(), -1),
    wallPatchOrNot_(mesh_.boundary().size(), 0),
    totalNCoarseFaces_(0),
    nLocalCoarseFaces_(0),
    nLocalWallCoarseFaces_(0),
    nLocalFineFaces_(0),
    totalNFineFaces_(0),
    constAlbedo_(false),
    timestepsInADay_(24),
    iterCounter_(0),
    pivotIndices_(0),
    FmyProc_(),
    globalFaceFacesLocal_(),
    maxJacobiIter_(200),
    jacobiTol_(1.0e-6)
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solarLoad::directAndDiffuse::~directAndDiffuse()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::solarLoad::directAndDiffuse::read()
{
    if (solarLoadModel::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


void Foam::solarLoad::directAndDiffuse::insertMatrixElements
(
    const globalIndex& globalNumbering,
    const label procI,
    const labelListList& globalFaceFaces,
    const scalarListList& viewFactors,
    scalarSquareMatrix& Fmatrix
)
{
    forAll(viewFactors, faceI)
    {
        const scalarList& vf = viewFactors[faceI];
        const labelList& globalFaces = globalFaceFaces[faceI];

        label globalI = globalNumbering.toGlobal(procI, faceI);
        forAll(globalFaces, i)
        {
            Fmatrix[globalI][globalFaces[i]] = vf[i];
        }
    }

    //Info << "Fmatrix: " << Fmatrix << endl;
}

void Foam::solarLoad::directAndDiffuse::insertScalarListListElements
(
    const globalIndex& globalNumbering,
    const label procI,
    const labelListList& sunskyMap,
    const labelListList& globalFaceFaces,
    const scalarListList& localCoeffs,
    scalarListList& globalCoeffs,
    const word& coarseOrFine
)
{
    forAll(localCoeffs, vectorId)
    {
        const scalarList& vf = localCoeffs[vectorId];

        forAll(vf, faceI)
        {
            if (coarseOrFine == "coarse")
            {
                globalCoeffs[vectorId][sunskyMap[procI][faceI]] = vf[faceI];
            }
            else
            {
                label globalI = globalNumbering.toGlobal(procI, faceI);
                globalCoeffs[vectorId][globalI] = vf[faceI];
            }
        }
    }
}

void Foam::solarLoad::directAndDiffuse::calculate()
{
    // Store previous iteration
    qs_.storePrevIter();

    scalarField compactCoarseA(map_->constructSize(), 0.0);
    scalarField compactCoarseHo(map_->constructSize(), 0.0);

    globalIndex globalNumbering(nLocalCoarseFaces_);

    // Fill local averaged Albedo(A) and external heatFlux(Ho)
    DynamicList<scalar> localCoarseAave(nLocalCoarseFaces_);
    DynamicList<scalar> localCoarseHoave(nLocalCoarseFaces_);

    volScalarField::Boundary& qsBf = qs_.boundaryFieldRef();

    forAll(selectedPatches_, i)
    {
        label patchID = selectedPatches_[i];

        const scalarField& sf = mesh_.magSf().boundaryField()[patchID];

        fvPatchScalarField& qsPatch = qsBf[patchID];

        solarLoadViewFactorFixedValueFvPatchScalarField& qsp =
            refCast
            <
                solarLoadViewFactorFixedValueFvPatchScalarField
            >(qsPatch);

        const scalarList ab = qsp.albedo();

        const scalarList& Hoi = qsp.qso();

        const polyPatch& pp = coarseMesh_.boundaryMesh()[patchID];
        const labelList& coarsePatchFace = coarseMesh_.patchFaceMap()[patchID];

        scalarList Aave(pp.size(), 0.0);
        scalarList Hoiave(Aave.size(), 0.0);

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
                // albedo and external flux area weighting
                forAll(fineFaces, j)
                {
                    label faceI = fineFaces[j];
                    Aave[coarseI] += (ab[faceI]*sf[faceI])/area;
                    Hoiave[coarseI] += (Hoi[faceI]*sf[faceI])/area;
                }
            }
        }

        //localCoarseTave.append(Tave);
        localCoarseAave.append(Aave);
        localCoarseHoave.append(Hoiave);
    }

    // Fill the local values to distribute
    SubList<scalar>(compactCoarseA,nLocalCoarseFaces_) = localCoarseAave;
    SubList<scalar>(compactCoarseHo,nLocalCoarseFaces_) = localCoarseHoave;

    // Distribute data
    map_->distribute(compactCoarseA);
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
    scalarField A(totalNCoarseFaces_, 0.0);
    scalarField qsExt(totalNCoarseFaces_, 0.0);

    // Fill lists from compact to global indexes.
    forAll(compactCoarseA, i)
    {
        A[compactGlobalIds[i]] = compactCoarseA[i];
        qsExt[compactGlobalIds[i]] = compactCoarseHo[i];
    }

    Pstream::listCombineGather(A, maxEqOp<scalar>());
    Pstream::listCombineGather(qsExt, maxEqOp<scalar>());

    Pstream::listCombineScatter(A);
    Pstream::listCombineScatter(qsExt);

    Time& time = const_cast<Time&>(mesh_.time());
    // Read sunPosVector list
    dictionary sunPosVectorIO;
    sunPosVectorIO.add(
        "file",
        fileName
        (
            mesh_.time().constant()
            /"sunPosVector"
        )
    );
    Function1s::TableFile<vector> sunPosVector
    (
        "sunPosVector",
        sunPosVectorIO
    );
    // look for the correct range
    label lo = 0;
    label hi = 0;
    scalarField sunPosVector_x = sunPosVector.x();
    forAll(sunPosVector_x, i)
    {
        if (time.value() >= sunPosVector_x[i])
        {
            lo = hi = i;
        }
        else
        {
            hi = i;
            break;
        }
    }
    scalar hi_fraction = 0;
    if (lo != hi) //if timestep is between two time values in sunPosVector
    {
        hi_fraction = (time.value() - sunPosVector_x[lo]) / (sunPosVector_x[hi] - sunPosVector_x[lo]);
    }

    // -------------------------------------------------------------------------
    // Parallel Jacobi solver for: C*q = b
    //   C(i,j) = delta_ij/(1-A_j) - (A_j/(1-A_j))*F(i,j)
    //   b(i)   = Isol(i) - sum(qsExt)   [derived by expanding the original
    //            O(N^2) loop: diagonal adds Isol[i]-qsExt[i], each off-diagonal
    //            adds -qsExt[j], total = Isol[i] - sum_j(qsExt[j])]
    //
    // Each proc owns a contiguous chunk [globalStart, globalStart+nLocal).
    // FmyProc_[localI] / globalFaceFacesLocal_[localI] hold the sparse row.
    // After each Jacobi step, q_new is allGathered so all procs hold the full
    // solution vector q for the next iteration's SpMV.
    // -------------------------------------------------------------------------
    const label myProc = Pstream::myProcNo();
    const label globalStart = globalNumbering.toGlobal(myProc, 0);
    const scalar sumQsExt = sum(qsExt);

    // --- Assemble local RHS (b) and diagonal C ---
    scalarField b_local(nLocalCoarseFaces_, 0.0);
    scalarField Cdiag_local(nLocalCoarseFaces_, 0.0);
    for (label localI = 0; localI < nLocalCoarseFaces_; localI++)
    {
        const label globalI = globalStart + localI;

        // Interpolated solar irradiance for this face and timestep
        scalar Isol =   skyViewCoeffGlobal_()[lo][localI]*(1 - hi_fraction)
                      + skyViewCoeffGlobal_()[hi][localI]*hi_fraction
                      + sunViewCoeffGlobal_()[lo][localI]*(1 - hi_fraction)
                      + sunViewCoeffGlobal_()[hi][localI]*hi_fraction;
        b_local[localI] = Isol - sumQsExt;

        // Diagonal of C: 1/(1-A_i) - A_i/(1-A_i)*F(i,i)
        // F(i,i) is 0 for convex/flat surfaces; search just in case.
        scalar Fdiag = 0.0;
        const labelList& cols = globalFaceFacesLocal_[localI];
        forAll(cols, k)
        {
            if (cols[k] == globalI) { Fdiag = FmyProc_[localI][k]; break; }
        }
        const scalar oneMinusAi = max(1.0 - A[globalI], VSMALL);
        Cdiag_local[localI] = 1.0/oneMinusAi - (A[globalI]/oneMinusAi)*Fdiag;
    }

    // --- Jacobi iteration ---
    scalarField q(totalNCoarseFaces_, 0.0);
    scalarField q_new(totalNCoarseFaces_, 0.0);

    scalar residual = GREAT;
    label iter = 0;
    while (iter < maxJacobiIter_ && residual > jacobiTol_)
    {
        // Each proc updates its own chunk of q_new; other entries stay 0.
        for (label localI = 0; localI < nLocalCoarseFaces_; localI++)
        {
            const label globalI = globalStart + localI;
            scalar rhs = b_local[localI];
            const labelList& cols = globalFaceFacesLocal_[localI];
            const scalarList& Frow = FmyProc_[localI];
            forAll(cols, k)
            {
                const label globalJ = cols[k];
                if (globalJ != globalI)
                {
                    const scalar oneMinusAj = max(1.0 - A[globalJ], VSMALL);
                    rhs += (A[globalJ]/oneMinusAj) * Frow[k] * q[globalJ];
                }
            }
            q_new[globalI] = rhs / Cdiag_local[localI];
        }

        // AllGather: each proc contributed only its chunk; sum = allgather
        Pstream::listCombineGather(q_new, sumOp<scalar>());
        Pstream::listCombineScatter(q_new);

        // Convergence check on local chunk (max absolute change)
        scalar localRes = 0.0;
        for (label localI = 0; localI < nLocalCoarseFaces_; localI++)
        {
            localRes = max(localRes, mag(q_new[globalStart + localI]
                                       - q[globalStart + localI]));
        }
        reduce(localRes, maxOp<scalar>());
        residual = localRes;

        Swap(q, q_new);
        q_new = scalar(0.0);
        iter++;
    }

    Info<< "\nJacobi solver: " << iter << " iterations, residual = "
        << residual << endl;



    label globCoarseId = 0;
    //label globFineId = 0;
    label fineFaceNo = 0;
    forAll(selectedPatches_, i)
    {
        const label patchID = selectedPatches_[i];
        const polyPatch& pp = mesh_.boundaryMesh()[patchID];
        if (pp.size() > 0)
        {
            scalarField& qsp = qsBf[patchID];
            const scalarField& sf = mesh_.magSf().boundaryField()[patchID];
            const labelList& agglom = finalAgglom_[patchID];
            label nAgglom = max(agglom)+1;

            labelListList coarseToFine(invertOneToMany(nAgglom, agglom));

            const labelList& coarsePatchFace =
                coarseMesh_.patchFaceMap()[patchID];

            scalar heatFlux = 0.0;
            forAll(coarseToFine, coarseI)
            {
                label globalCoarse =
                    globalNumbering.toGlobal(Pstream::myProcNo(), globCoarseId);
                const label coarseFaceID = coarsePatchFace[coarseI];
                const labelList& fineFaces = coarseToFine[coarseFaceID];
                forAll(fineFaces, k)
                {
                    label faceI = fineFaces[k];

                    qsp[faceI] = q[globalCoarse];
                    if (isA<wallFvPatch>(mesh_.boundary()[patchID]))
                    {
                        qsp[faceI] -= (sunViewCoeffGlobal_()[lo][globCoarseId]*(1-hi_fraction) + sunViewCoeffGlobal_()[hi][globCoarseId]*(hi_fraction)) * (1-A[globalCoarse]);
                        qsp[faceI] += (solarLoadFineFacesGlobal_()[lo][fineFaceNo+faceI]*(1-hi_fraction) + solarLoadFineFacesGlobal_()[hi][fineFaceNo+faceI]*(hi_fraction)) * (1-A[globalCoarse]);
                    }
                    heatFlux += qsp[faceI]*sf[faceI];
                }
                globCoarseId ++;
            }
        }
        if (isA<wallFvPatch>(mesh_.boundary()[patchID]))
        {
            fineFaceNo += pp.size();
        }
    }

    if (debug)
    {
        forAll(qsBf, patchID)
        {
            const scalarField& qsp = qs_.boundaryField()[patchID];
            const scalarField& magSf = mesh_.magSf().boundaryField()[patchID];
            const scalar heatFlux = gSum(qsp*magSf);
            Info<< "Total heat transfer rate at patch: "
                << patchID << " "
                << heatFlux << endl;
        }
    }

    // Relax qs if necessary
    qs_.relax();
}


Foam::tmp<Foam::volScalarField> Foam::solarLoad::directAndDiffuse::Rp() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Rp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar
            (
                "zero",
                dimMass/pow3(dimTime)/dimLength/pow4(dimTemperature),
                0.0
            )
        )
    );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::solarLoad::directAndDiffuse::Ru() const
{
    return tmp<volScalarField::Internal>
    (
        new volScalarField::Internal
        (
            IOobject
            (
                "Ru",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimMass/dimLength/pow3(dimTime), 0.0)
        )
    );
}

// ************************************************************************* //
