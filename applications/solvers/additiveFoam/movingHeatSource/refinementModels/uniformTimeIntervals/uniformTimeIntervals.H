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
    Foam::refinementModels::uniformTimeIntervals

Description
    This refinement control scheme divides the total scan-path duration into
    a fixed number of uniform time intervals. At each interval, the refinement
    marker field is updated using the current temperature field and the
    projected scan path over the next interval.

SourceFiles
    uniformTimeIntervals.C

\*---------------------------------------------------------------------------*/

#ifndef uniformTimeIntervals_H
#define uniformTimeIntervals_H

#include "refinementModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace refinementModels
{

/*---------------------------------------------------------------------------*\
                       Class uniformTimeIntervals
\*---------------------------------------------------------------------------*/

class uniformTimeIntervals
:
    public Foam::refinementModel
{
    // Private Data

        //- Coefficients dictionary
        dictionary coeffs_;


protected:

    // Protected Data

        //- Number of uniform scan-path intervals
        label nIntervals_;

        //- Duration of each interval
        scalar intervalLength_;

        //- Time at which to update refinement again
        scalar updateTime_;


public:

    //- Runtime type information
    TypeName("uniformTimeIntervals");


    // Constructors

        //- Construct from components
        uniformTimeIntervals
        (
            const PtrList<heatSourceModel>& sources,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~uniformTimeIntervals()
    {}


    // Member Functions

        //- Update the refinement marker field
        virtual bool update();

        //- Read the heatSourceProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace refinementModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
