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
    Foam::heatSourceModels::modifiedSuperGaussian

Description
    modifiedSuperGaussian model.

    Reference:
    \verbatim
   \endverbatim

SourceFiles
    modifiedSuperGaussian.C

\*---------------------------------------------------------------------------*/

#ifndef modifiedSuperGaussian_H
#define modifiedSuperGaussian_H

#include "heatSourceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heatSourceModels
{

/*---------------------------------------------------------------------------*\
                       Class modifiedSuperGaussian
\*---------------------------------------------------------------------------*/

class modifiedSuperGaussian
:
    public heatSourceModel
{
    // Private Data
    
        //- Pointer to mesh information
        const fvMesh& mesh_;
        
        //- Distribution shape coefficient of modifiedSuperGaussian distribution
        scalar k_;

        //- Geometric shape coefficient of modifiedSuperGaussian distribution
        scalar m_;

public:

    //- Runtime type information
    TypeName("modifiedSuperGaussian");


    // Constructors

        //- Construct from components
        modifiedSuperGaussian
        (
            const word& sourceName,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~modifiedSuperGaussian()
    {}


    // Member Functions

        //- Return the source volumetric heating
        virtual tmp<volScalarField> qDot();
        
        //- Return the D4sigma value of the modifiedSuperGaussian
        virtual scalar D4sigma() const;

        //- Read the heatSourceProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heatSourceModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
