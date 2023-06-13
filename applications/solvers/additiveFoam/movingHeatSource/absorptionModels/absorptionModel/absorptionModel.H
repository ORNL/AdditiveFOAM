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
    Foam::absorptionModel

Description
    Abstract base class for absorption models

SourceFiles
    absorptionModel.C

\*---------------------------------------------------------------------------*/

#ifndef absorptionModel_H
#define absorptionModel_H

#include "dictionary.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class absorptionModel Declaration
\*---------------------------------------------------------------------------*/

class absorptionModel
:
    public IOdictionary
{
    // Private member functions
    
        //- Construct the base IO object
        IOobject createIOobject
        (
            const word& sourceName,
            const dictionary& dict,
            const fvMesh& mesh
        ) const;
    
protected:

    // Protected data

        //- Name of beam
        const word& sourceName_;

        //- Main heat source dictionary
        const dictionary& heatSourceDict_;
           
           //- Pointer to mesh info
        const fvMesh& mesh_;

        //- Individual beam dictionary
        dictionary beamDict_;

        //- Absorption model coefficient dictionary
        dictionary absorptionModelCoeffs_;
        

public:

    //- Runtime type information
    TypeName("absorptionModel");
    
    //- Default heatSourceProperties dictionary name
    static const word heatSourceDictName;
    
    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        absorptionModel,
        dictionary,
        (
            const word& sourceName,
            const dictionary& dict,
            const fvMesh& mesh
        ),
        (sourceName, dict, mesh)
    );

    // Constructors

        //- Construct
        absorptionModel
        (
            const word& type,
            const word& sourceName,
            const dictionary& dict,
            const fvMesh& mesh
        );
        
        //- Disallow default bitwise copy construction
        absorptionModel(const absorptionModel&);
        
    // Selectors

        //- Return a reference to the selected absorption model
        static autoPtr<absorptionModel> New
        (
            const word& sourceName,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~absorptionModel()
    {}


    // Member Functions

        //- Return const-access to the absorption
        virtual scalar eta
        (
            const scalar& aspectRatio
        ) const = 0;

        //- Read the absorption properties
        virtual bool read();
        
    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const absorptionModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
