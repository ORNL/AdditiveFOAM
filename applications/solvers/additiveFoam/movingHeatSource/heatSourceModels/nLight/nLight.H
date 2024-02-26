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
    Foam::heatSourceModels::nLight

Description
    Beam profile based on nLight system described in the listed reference.
    The article describes a 2D profile which was extended to a generalized 3D
    distribution by Jamie Stump and John Coleman for use with additiveFoam.
    The exact profiles given in the nLight article may be approximated by
    fitting the model coefficients in this class.

    Reference:
        Vanmunster, L., Goossens, L., Kinds, Y., van Hooreweder, B.,
        & Vrancken, B. (2023). Effect of a variable laser beam profile 
        system on productivity and surface quality of 316L stainless steel 
        parts produced by Laser Powder Bed Fusion. CIRP Annals.
        https://doi.org/10.1016/j.cirp.2023.04.035


SourceFiles
    nLight.C

\*---------------------------------------------------------------------------*/

#ifndef nLight_H
#define nLight_H

#include "heatSourceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heatSourceModels
{

/*---------------------------------------------------------------------------*\
                                  Class nLight
\*---------------------------------------------------------------------------*/

class nLight
:
    public heatSourceModel
{
    // Private Data
    
        //- Pointer to mesh information
        const fvMesh& mesh_;
        
        //- Spot/beam weight
        scalar alpha_;
        
        //- Dictionary for spot shape parameters
        const dictionary& spotCoeffs_;
    
        //- Laser mode [0-6]
        scalar mode_;
    
        //- Spot distribution shape coefficient
        scalar ks_;

        //- Spot geometric shape coefficient
        scalar ms_;
        
        //- Spot dimensions
        vector ds_;
        
        //- Spot area
        scalar As_;
        
        //- Dictionary for ring shape parameters
        const dictionary& ringCoeffs_;
        
        //- Ring distribution shape coefficient
        scalar kr_;
        
        //- Ring geometric shape coefficient
        scalar mr_;
        
        //- Ring center radius
        scalar R_;
        
        //- Ring half width
        scalar r_;
        
        //- Ring area
        scalar Ar_;


public:

    //- Runtime type information
    TypeName("nLight");


    // Constructors

        //- Construct from components
        nLight
        (
            const word& sourceName,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~nLight()
    {}


    // Member Functions

        inline virtual scalar weight(const vector& d);

        inline virtual dimensionedScalar V0();
        
        // Normalization factors for spot and ring integration
        // NOTE: THESE DO NOT INCLUDE THE BEAM DEPTH!
        inline virtual dimensionedScalar As();
        inline virtual dimensionedScalar Ar();

        //- Read the heatSourceProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heatSourceModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
