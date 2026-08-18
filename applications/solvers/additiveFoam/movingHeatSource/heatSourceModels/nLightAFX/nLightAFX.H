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
    Foam::heatSourceModels::nLightAFX

Description
    Two-component nLight AFX heat source model with a shared projected axial
    distribution applied in the negative z direction from the beam plane. The
    beam plane is expected to coincide with the material surface.

SourceFiles
    nLightAFX.C

\*---------------------------------------------------------------------------*/

#ifndef nLightAFX_H
#define nLightAFX_H

#include "heatSourceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heatSourceModels
{

/*---------------------------------------------------------------------------*\
                          Class nLightAFX
\*---------------------------------------------------------------------------*/

class nLightAFX
:
    public heatSourceModel
{
    // Private Data

        //- Pointer to mesh information
        const fvMesh& mesh_;

        //- Fraction of power applied to Gaussian component 1
        scalar alpha_;

        //- Radius of Gaussian component 0
        scalar r0_;

        //- Radius of Gaussian component 1
        scalar r1_;

        //- Standard deviation of Gaussian component 0
        scalar sigma0_;

        //- Standard deviation of Gaussian component 1
        scalar sigma1_;

        //- nSlope in n = nSlope*log2(x) + nIntercept
        scalar nSlope_;

        //- nIntercept in n = nSlope*log2(x) + nIntercept
        scalar nIntercept_;

        //- Axial decay exponent
        scalar k_;

        //- One-sided volume integral of Gaussian component 0
        scalar a0_;

        //- One-sided volume integral of Gaussian component 1
        scalar a1_;

public:

    //- Runtime type information
    TypeName("nLightAFX");


    // Constructors

        //- Construct from components
        nLightAFX
        (
            const word& sourceName,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~nLightAFX()
    {}


    // Member Functions

        inline virtual scalar weight(const vector& d);

        //- One-sided volume integral of a projected Gaussian ring
        inline scalar ai(scalar x, scalar s, scalar d, scalar k);

        inline virtual dimensionedScalar V0();

        //- Read the heatSourceDict dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heatSourceModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
