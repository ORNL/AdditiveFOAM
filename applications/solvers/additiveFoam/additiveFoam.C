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

Application
    additiveFoam

Description
    A transient heat transfer and fluid flow solver for additive manufacturing
    simulations.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "zeroGradientFvPatchFields.H"
#include "IFstream.H"
#include "uniformDimensionedFields.H"
#include "pressureReference.H"
#include "findRefCell.H"

#include "fvmDiv.H"
#include "fvmDdt.H"
#include "fvcSurfaceIntegrate.H"
#include "fvcVolumeIntegrate.H"
#include "fvmLaplacian.H"
#include "constrainPressure.H"
#include "constrainHbyA.H"
#include "adjustPhi.H"
#include "pimpleControl.H"
#include "fvCorrectPhi.H"
#include "Polynomial.H"

#include "EulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"

#include "additiveFoamInfo.H"
#include "movingHeatSourceModel.H"
#include "thermoPath.H"
#include "interpolateXY.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    using namespace Foam;

    #include "postProcess.H"
    #include "setRootCase.H"

    Foam::AdditiveFoamInfo::write();

    #include "createTime.H"
    #include "createMesh.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Initialize time-stepping controls
    scalar DiNum = 0.0;
    scalar alphaCoNum = 0.0;
    movingHeatSourceModel sources(mesh);

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "updateProperties.H"

        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        sources.update();

        mesh.update();

        runTime++;

        Info<< "Time = " << runTime.name() << nl << endl;

        #include "solutionControls.H"

        while (pimple.loop())
        {
            #include "moveMesh.H"

            if (fluidInDomain)
            {
                #include "pU/UEqn.H"
                #include "pU/pEqn.H"
            }
        }

        #include "thermo/TEqn.H"

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    return 0;
}

// ************************************************************************* //
