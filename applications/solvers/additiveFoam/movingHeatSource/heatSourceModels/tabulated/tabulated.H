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
    Tabulated 2D heat source model with a projected axial distribution applied
    in the negative z direction from the beam plane. The beam plane is expected
    to coincide with the material surface.

SourceFiles
    tabulated.C

\*---------------------------------------------------------------------------*/

#ifndef tabulated_H
#define tabulated_H

#include "heatSourceModel.H"
#include "tabulatedProfile.H"

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

        //- nSlope in n = nSlope*log2(a) + nIntercept
        scalar nSlope_;

        //- nIntercept in n = nSlope*log2(a) + nIntercept
        scalar nIntercept_;

        //- Axial decay exponent
        scalar k_;

        //- Planar profile and derived beam metrics
        tabulatedProfile profile_;


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

        virtual void update();

        inline virtual scalar weight(const vector& r) const;

        //- Read the heatSourceDict dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heatSourceModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
