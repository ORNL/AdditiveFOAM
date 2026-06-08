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
    nLight AFX heat source model for additive manufacturing

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

        //- Fraction of power applied to the outer component
        scalar alpha_;

        //- Radius of the inner Gaussian ring
        scalar r0_;

        //- Radius of the outer Gaussian ring
        scalar r1_;

        //- Standard deviation of the inner Gaussian ring
        scalar sigma0_;

        //- Standard deviation of the outer Gaussian ring
        scalar sigma1_;

        //- Inner axial decay coefficient in the form n = Alog2(x) + B
        scalar A0_;

        //- Inner axial decay coefficient in the form n = Alog2(x) + B
        scalar B0_;

        //- Outer axial decay coefficient in the form n = Alog2(x) + B
        scalar A1_;

        //- Outer axial decay coefficient in the form n = Alog2(x) + B
        scalar B1_;

        //- Inner axial decay exponent
        scalar n0_;

        //- Outer axial decay exponent
        scalar n1_;

        //- Inner normalization factor
        scalar a0_;

        //- Outer normalization factor
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

        inline scalar ai(scalar x, scalar s, scalar d, scalar n);

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
