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

\*---------------------------------------------------------------------------*/

#include "uniformTimeIntervals.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace refinementModels
{
    defineTypeNameAndDebug(uniformTimeIntervals, 0);

    addToRunTimeSelectionTable
    (
        refinementModel,
        uniformTimeIntervals,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementModels::uniformTimeIntervals::uniformTimeIntervals
(
    const PtrList<heatSourceModel>& sources,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    Foam::refinementModel(typeName, sources, dict, mesh),
    coeffs_(refinementDict_.optionalSubDict(typeName + "Coeffs")),
    nIntervals_(coeffs_.lookup<label>("intervals")),
    intervalLength_(Zero),
    updateTime_(Zero)
{
    //- Set interval length from the active scan-path duration
    intervalLength_ =
        max(scalar(0), scanEndTime_ - mesh.time().startTime().value())
       /nIntervals_;

    Info<< typeName << ": Interval duration set to "
        << intervalLength_ << " s." << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::refinementModels::uniformTimeIntervals::update()
{
    if ((updateTime_ - mesh_.time().value()) < small)
    {
        //- Update next refinement time
        updateTime_ = mesh_.time().value() + intervalLength_;

        //- Reactive refinement based on temperature
        Foam::refinementModel::markTemperature();

        //- Scan path completed
        if ((scanEndTime_ - mesh_.time().value()) < small)
        {
            Info<< typeName << ": "
                << "Scan path completed. "
                << "Continuing AMR checks for possible mesh coarsening."
                << endl;

            updateTime_ = mesh_.time().value() + intervalLength_;

            return true;
        }

        Info<< typeName << ": Updating AMR marker field." << endl;

        //- Predictive refinement along scan path over the next interval
        Foam::refinementModel::markScanPathTime(updateTime_);
    }

    return true;
}


bool Foam::refinementModels::uniformTimeIntervals::read()
{
    if (Foam::refinementModel::read())
    {
        coeffs_ = refinementDict_.optionalSubDict(typeName + "Coeffs");

        coeffs_.lookup("intervals") >> nIntervals_;

        intervalLength_ =
            max(scalar(0), scanEndTime_ - mesh_.time().startTime().value())
           /nIntervals_;

        return true;
    }

    return false;
}


// ************************************************************************* //
