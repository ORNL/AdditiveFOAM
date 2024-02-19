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
    Foam::refinementControllers::refinementController

Description
    Abstract base class for refinementController models

SourceFiles
    refinementController.C

\*---------------------------------------------------------------------------*/

#ifndef refinementController_H
#define refinementController_H

#include "labelVector.H"
#include "hexMatcher.H"
#include "dictionary.H"
#include "volFields.H"
#include "movingBeam.H"
#include "treeBoundBox.H"
#include "heatSourceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class refinementController Declaration
\*---------------------------------------------------------------------------*/

class refinementController
:
    public IOdictionary
{
    // Private member functions
    
        //- Construct the base IO object
        IOobject createIOobject
        (
            const dictionary& dict,
            const fvMesh& mesh
        ) const;
    
protected:

    // Protected data
    
        //- PtrList to heatSourceModels
        const PtrList<heatSourceModel>& sources_;

        //- Pointer to mesh information
        const fvMesh& mesh_;

        //- Main heat source dictionary
        const dictionary& heatSourceDict_;
                
        //- Sub-dictionary for refinementController entries
        dictionary refinementDict_;
        
        //- Switch for refinementField generation
        const bool refine_;
                
        //- Number of refinement levels
        label nLevels_;

        //- Time index of last refinement operation
        label lastRefinementIndex_;

        //- Marker field for refinement operations
        volScalarField refinementField_;

public:

    //- Runtime type information
    TypeName("refinementController");
    
    //- Default heatSourceProperties dictionary name
    static const word refinementControllerDictName;
    
    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        refinementController,
        dictionary,
        (
            const PtrList<heatSourceModel>& sources,
            const dictionary& dict,
            const fvMesh& mesh
        ),
        (sources, dict, mesh)
    );


    // Constructors

        //- Construct for mixture
        refinementController
        (
            const word& type,
            const PtrList<heatSourceModel>& sources,
            const dictionary& dict,
            const fvMesh& mesh
        );
        
        //- Disallow default bitwise copy construction
        refinementController(const refinementController&);
        
    // Selectors

        //- Return a reference to the selected heatSource model
        static autoPtr<refinementController> New
        (
            const PtrList<heatSourceModel>& sources,
            const dictionary& dict,
            const fvMesh& mesh
        );

    //- Destructor
    virtual ~refinementController()
    {}


    // Member Functions
    
        //- Return the refinement marker field
        virtual volScalarField refinementField()
        {
            return refinementField_;
        };
        
        //- Set refinementField using temperature constraints
        virtual void setRefinementField();
    
        //- Update refinementField and return true at specified time
        virtual bool update() = 0;

        //- Read the heat source dictionary
        virtual bool read() = 0;
        
        
    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const refinementController&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
