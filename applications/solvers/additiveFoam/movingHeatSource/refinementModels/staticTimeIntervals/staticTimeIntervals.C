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

\*---------------------------------------------------------------------------*/

#include "staticTimeIntervals.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace refinementModel
{
    defineTypeNameAndDebug(staticTimeIntervals, 0);
    addToRunTimeSelectionTable
    (
        refinementModel,
        staticTimeIntervals,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementModel::staticTimeIntervals::staticTimeIntervals
(
    const PtrList<heatSourceModel>& sources,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    refinementModel(typeName, sources, dict, mesh),
    coeffs_(refinementDict_.optionalSubDict(typeName + "Coeffs")),
    intervals_(coeffs_.lookup<scalar>("intervals")),
    intervalLength_(0.0),
    updateTime_(0.0)
{
    //- Set interval time to the number of intervals between start and end
    intervalLength_ =
        max(0, endTime_ - mesh.time().startTime().value()) / intervals_;

    Info << "staticTimeIntervals: Interval time duration set to " 
         << intervalLength_ << " s." << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::refinementModel::staticTimeIntervals::update()
{
    if (updateTime_ - mesh_.time().value() < small)
    {
        // Update next refinement time
        updateTime_ = mesh_.time().value() + intervalLength_;

        //- Refine in regions above specified temperature
        refinementModel::refineUsingTemperature();

        //- Don't perform additional refinements if scan path is completed
        if ((endTime_ - mesh_.time().value()) < small)
        {
            Info << "staticTimeIntervals: Scan path completed. Continuing AMR"
                 << " checks for possible mesh coarsening" << endl;
                 
            updateTime_ = mesh_.time().value() + intervalLength_;
                 
            return true;
        }
        
        Info << "staticTimeIntervals: Updating AMR marker field." << endl;
        
        refinementModel::refineUsingTime(updateTime_);
    }

    return true;
}


bool Foam::refinementModel::staticTimeIntervals::read()
{
    if (refinementModel::read())
    {
        refinementDict_ = optionalSubDict(type() + "Coeffs");

        //- Mandatory entries
        refinementDict_.lookup("intervals") >> intervals_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
