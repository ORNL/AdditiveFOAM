/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::mixedTemperatureFvPatchScalarField

Description
    This boundary condition provides a mixedTemperature condition,
    calculated as:

        \f[
            Q = h * (T - T_inf) + e * sigma * (T^4 - T_inf^$)
        \f]

    where
    \vartable
        h      | Convection coefficient [W/m^2/K]
        e      | Surface emissivity [-]
        sigma  | Stefan-Boltzmann constant [5.67e-8 W/m^2/K^4]
    \endtable

Usage
    \table
        Property   | Description                  | Req'd? | Default
        h          | Uniform convection coeff     | yes    |
        emissivity | Uniform surface emissivity   | no     | 0.0
        Tinf       | Ambient temperature          | yes    |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type       mixedTemperature;
        h          10.0;
        emissivity 0.4;
    }
    \endverbatim

SourceFiles
    mixedTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef BC_H
#define BC_H

#include "mixedFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class mixedTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class mixedTemperatureFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private Data

        //- Convective heat transfer coefficient [W/m^2/K]
        scalar h_;
        
        //- Effective emissivity of the boundary
        scalar emissivity_;
        
        //- Ambient temperature
        scalarField Tinf_;


public:

    //- Runtime type information
    TypeName("mixedTemperature");


    // Constructors

        //- Construct from patch, internal field and dictionary
        mixedTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given mixedTypeFvPatchField
        //  onto a new patch
        mixedTemperatureFvPatchScalarField
        (
            const mixedTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        mixedTemperatureFvPatchScalarField
        (
            const mixedTemperatureFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        mixedTemperatureFvPatchScalarField
        (
            const mixedTemperatureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new mixedTemperatureFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        // Mapping functions

            //- Map the given fvPatchField onto this fvPatchField
            virtual void map
            (
                const fvPatchScalarField&,
                const fvPatchFieldMapper&
            );

            //- Reset the fvPatchField to the given fvPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const fvPatchScalarField&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
