/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                Copyright (C) 2023 Oak Ridge National Laboratory                
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
    Foam::absorptionModel::Kelly

Description
    Kelly absorption model for the optical absorption due to multiple internal
    reflections in simple geometric shapes (e.g. cone and cylinder).

    Reference:
    \verbatim
        F. J. Kelly,
        "An Equation for the Local Thermal Emissivity at the
        Vertex of a Diffuse Conical or V-Groove Cavity",
        Applied Optics, Vol. 5, No. 6, June 1966.
   \endverbatim

SourceFiles
    KellyAbsorption.C

\*---------------------------------------------------------------------------*/

#ifndef KellyAbsorption_H
#define KellyAbsorption_H

#include "absorptionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace absorptionModels
{

/*---------------------------------------------------------------------------*\
                       Class Kelly
\*---------------------------------------------------------------------------*/

class Kelly
:
    public absorptionModel
{
    // Private Data
    
        //- Conical or cylindrical
        word geometry_;

        //- Constant absorption scalar
        dimensionedScalar eta0_;
        
        //- Minimum absorption
        dimensionedScalar etaMin_;

public:

    //- Runtime type information
    TypeName("Kelly");


    // Constructors

        //- Construct from components
        Kelly
        (
            const word& sourceName,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~Kelly()
    {}

    // Member Functions

        //- Return the absorption as a scalar field
        virtual scalar eta
        (
            const scalar& aspectRatio
        ) const;

        //- Read the absorption properties
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace absorptionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
