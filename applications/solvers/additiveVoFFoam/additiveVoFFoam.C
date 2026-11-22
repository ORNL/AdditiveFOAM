/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                Copyright (C) 2026 Oak Ridge National Laboratory
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
    additiveVoFFoam

Description
    A transient heat transfer and fluid flow solver for additive manufacturing
    simulations with a metal-gas free surface. The metal is captured with a
    volume of fluid method and the beam is absorbed at the first metal surface
    it reaches.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "zeroGradientFvPatchFields.H"
#include "uniformDimensionedFields.H"
#include "pressureReference.H"
#include "findRefCell.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "fvcAverage.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcFlux.H"
#include "fvcGrad.H"
#include "fvcMeshPhi.H"
#include "fvcReconstruct.H"
#include "fvcSnGrad.H"
#include "fvcSurfaceIntegrate.H"
#include "constrainPressure.H"
#include "constrainHbyA.H"
#include "adjustPhi.H"
#include "pimpleControl.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "gaussConvectionScheme.H"
#include "upwind.H"
#include "mathematicalConstants.H"
#include "physicoChemicalConstants.H"

#include "additiveFoamInfo.H"
#include "movingHeatSources.H"
#include "solidLiquidThermo.H"
#include "thermo/sensibleEnergy.H"
#include "phaseChangeSolver.H"

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

    if (mesh.dynamic())
    {
        FatalErrorInFunction
            << "additiveVoFFoam requires a static mesh" << exit(FatalError);
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    movingHeatSources sources(mesh);

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        sources.update(alpha1);

        runTime++;

        Info<< "Time = " << runTime.name() << nl << endl;

        #include "solutionControls.H"

        while (pimple.loop())
        {
            #include "alphaEqn.H"

            #include "thermo/TEqn.H"
            #include "pU/UEqn.H"

            while (pimple.correct())
            {
                #include "pU/pEqn.H"
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
