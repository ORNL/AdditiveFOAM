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
    Foam::heatSourceModels::projectedGaussian

Description
    projectedGaussian heat source model for additive manufacting

SourceFiles
    projectedGaussian.C

\*---------------------------------------------------------------------------*/

#ifndef projectedGaussian_H
#define projectedGaussian_H

#include "heatSourceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heatSourceModels
{

/*---------------------------------------------------------------------------*\
                       Class projectedGaussian
\*---------------------------------------------------------------------------*/

class projectedGaussian
:
    public heatSourceModel
{
    // Private Data
    
        //- Pointer to mesh information
        const fvMesh& mesh_;
        
        //- weights for axial decay of 2D gaussian distribution in the form:
        //- n = Alog2(x) + B
        scalar A_;
        scalar B_;
        
        scalar k_;

public:

    //- Runtime type information
    TypeName("projectedGaussian");


    // Constructors

        //- Construct from components
        projectedGaussian
        (
            const word& sourceName,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~projectedGaussian()
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
