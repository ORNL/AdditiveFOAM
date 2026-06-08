/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                Copyright (C) 2023-2026 Oak Ridge National Laboratory
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
    Foam::heatSourceModels::tabulated

Description
    Tabulated 2D heat source model with a projected axial distribution

SourceFiles
    tabulated.C

\*---------------------------------------------------------------------------*/

#ifndef tabulated_H
#define tabulated_H

#include "heatSourceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heatSourceModels
{

/*---------------------------------------------------------------------------*\
                          Class tabulated
\*---------------------------------------------------------------------------*/

class tabulated
:
    public heatSourceModel
{
    // Private Data

        //- Pointer to mesh information
        const fvMesh& mesh_;

        //- Axial decay coefficient in the form n = Alog2(x) + B
        scalar A_;

        //- Axial decay coefficient in the form n = Alog2(x) + B
        scalar B_;

        //- Axial decay exponent
        scalar k_;

        //- Minimum x coordinate of the table
        scalar x0_;

        //- Minimum y coordinate of the table
        scalar y0_;

        //- Uniform x spacing
        scalar dx_;

        //- Uniform y spacing
        scalar dy_;

        //- Inverse uniform x spacing
        scalar invDx_;

        //- Inverse uniform y spacing
        scalar invDy_;

        //- Number of x points
        label nx_;

        //- Number of y points
        label ny_;

        //- Tabulated planar intensity, stored as i + nx*j
        scalarField f_;

        //- Integrated planar distribution
        scalar Axy_;


    // Private Member Functions

        //- Read tabulated planar distribution
        void readTable(const fileName& tableFile);

        //- Integrate tabulated planar distribution
        void integrateTable();


public:

    //- Runtime type information
    TypeName("tabulated");


    // Constructors

        //- Construct from components
        tabulated
        (
            const word& sourceName,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~tabulated()
    {}


    // Member Functions

        inline virtual scalar weight(const vector& d);

        inline virtual dimensionedScalar V0();

        //- Read the heatSourceProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heatSourceModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
