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

Class
    Foam::radiationModels::viewFactorSky

Description
    View factor radiation model.
    Sky temperature is used for non-wall patches, e.g. side and top boundaries.
    Written by Aytac Kubilay based on the original viewFactor model.
    
    The system solved is: C q = b
    where:
            Cij  = deltaij/Ej - (1/Ej - 1)Fij
            q    = heat flux
            b    = A eb - Ho
    and:
            eb   = sigma*T^4
            Ej   = emissivity
            Aij  = deltaij - Fij
            Fij  = view factor matrix


SourceFiles
    viewFactorSky.C

\*---------------------------------------------------------------------------*/

#ifndef radiationModelviewFactorSky_H
#define radiationModelviewFactorSky_H

#include "radiationModel.H"
#include "singleCellFvMesh.H"
#include "globalIndex.H"
#include "scalarListIOList.H"
#include "mapDistribute.H"
#include "volFields.H"
#include "hashedWordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{

/*---------------------------------------------------------------------------*\
                           Class viewFactorSky Declaration
\*---------------------------------------------------------------------------*/

class viewFactorSky
:
    public radiationModel
{
    // Private Data

        //- Agglomeration List
        labelListIOList finalAgglom_;

        //- Map distributed
        autoPtr<mapDistribute> map_;

        //- Coarse mesh
        singleCellFvMesh coarseMesh_;

        //- Net radiative heat flux [W/m^2]
        volScalarField qr_;

        //- Sparse view factor matrix (master only): per-row column indices and values
        labelListList FsparseCols_;
        scalarListList FsparseVals_;

        //- Selected patches
        labelList selectedPatches_;

        //- Total global coarse faces
        label totalNCoarseFaces_;

        //- Total local coarse faces
        label nLocalCoarseFaces_;

        //- Global index of the first local row (= globalNumbering.toGlobal(myProc,0))
        label globalOffset_;

        //- GMRES solver parameters
        label gmresRestart_;
        label gmresMaxIter_;
        scalar gmresTol_;

        //- Previous solution for warm-start across time steps
        scalarField qPrev_;

        //- List of grass patches
        hashedWordList grassPatches;


    // Private Member Functions

        //- Initialise
        void initialise();

        //- Sparse matrix-vector product: w = C*x
        //  C(i,j) = delta_ij/E_j - (1/E_j-1)*F_ij
        scalarField matvecC
        (
            const scalarField& x,
            const scalarField& E
        ) const;

        //- Restarted GMRES(gmresRestart_) solver for C*x = b
        //  x is in/out (warm start); returns final absolute residual norm
        scalar solveGMRES
        (
            const scalarField& b,
            scalarField& x,
            const scalarField& E
        ) const;


public:

    //- Runtime type information
    TypeName("viewFactorSky");


    // Constructors

        //- Construct from components
        viewFactorSky(const volScalarField& T);

        //- Construct from components
        viewFactorSky(const dictionary& dict, const volScalarField& T);

        //- Disallow default bitwise copy construction
        viewFactorSky(const viewFactorSky&) = delete;


    //- Destructor
    virtual ~viewFactorSky();


    // Member Functions

        // Edit

            //- Solve system of equation(s)
            void calculate();

            //- Read radiation properties dictionary
            bool read();

            //- Source term component (for power of T^4)
            virtual tmp<volScalarField> Rp() const;

            //- Source term component (constant)
            virtual tmp<volScalarField::Internal> Ru() const;


    // Access

        //- Const access to total radiative heat flux field
        inline const volScalarField& qr() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const viewFactorSky&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "viewFactorSkyI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
