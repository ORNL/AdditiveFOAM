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

Application
    additiveFoam

Description
    A transient heat transfer and fluid flow solver for additive manufacturing
    simulations.
    
\*---------------------------------------------------------------------------*/

#include "Make/additiveFoamVersion.H"

#include "fvCFD.H"
#include "pimpleControl.H"
#include "graph.H"
#include "Polynomial.H"
#include "interpolateXY/interpolateXY.H"
#include "movingHeatSourceModel.H"
#include "EulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{       
    #include "postProcess.H"    
    #include "setRootCase.H"
    
    Info << "AdditiveFOAM Information:" << nl;
    #ifdef ADDITIVEFOAM_VERSION
        Info << "Version:   " << ADDITIVEFOAM_VERSION << nl;
    #endif
    #ifdef ADDITIVEFOAM_GIT_DESCRIBE
        Info << "Build:     " << ADDITIVEFOAM_GIT_DESCRIBE << nl;
    #endif
    #ifdef ADDITIVEFOAM_GIT_SHA1
        Info << "Git SHA1:  " << ADDITIVEFOAM_GIT_SHA1 << nl << endl;
    #endif
    
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createTimeControls.H"
    #include "initContinuityErrs.H"
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    scalar DiNum = 0.0;
    scalar alphaCoNum = 0.0;
    movingHeatSourceModel sources(mesh);

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "updateProperties.H"

        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        sources.update();
        
        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "solutionControls.H"
        
        while (pimple.loop() && fluidInDomain)
        {
            #include "pU/UEqn.H"
            #include "pU/pEqn.H"
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
