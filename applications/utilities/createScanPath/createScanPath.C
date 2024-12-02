/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    createScanPath

\*---------------------------------------------------------------------------*/

#include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    #include "setRootCase.H"
    #include "createTime.H"
    
    IOdictionary dict
    (
        IOobject
        (
            "createScanPathDict",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    point2D p0 = dict.lookup<point2D>("minPoint");
    Point minPoint(p0.x(), p0.y());

    point2D p1 = dict.lookup<point2D>("maxPoint");
    Point maxPoint(p1.x(), p1.y());   
    
    scalar angle = dict.lookup<scalar>("angle");
    scalar hatch = dict.lookup<scalar>("hatch");
    
    label nRotations = dict.lookup<label>("nRotations");
    
    scalar power = dict.lookup<scalar>("power");
    scalar speed = dict.lookup<scalar>("speed");
    scalar dwellTime = dict.lookup<scalar>("dwellTime");
    
    bool biDirection = dict.lookupOrDefault<bool>("biDirection", true);
    
    // Create bounding box for scan vectors
    BoundBox boundingBox( minPoint, maxPoint );

    // Rotate and write scan vectors to file
    scalar rotation = 0.0;
    
    for (label i = 0; i < nRotations; ++i)
    {
        Info << "Creating scan path for rotation angle: " << rotation << endl;
        
        Path path( boundingBox, hatch, rotation );
        
        path.power = power;
        path.speed = speed;
        path.dwellTime = dwellTime;

        std::ostringstream oss;
        oss << std::fixed << std::setprecision(0) << rotation;
        std::string rotationString = oss.str();
        
        const fileName filename
        (
            runTime.constant()/"scanPath_" + std::to_string(i)
        );

        path.write(filename, biDirection);

        rotation += angle;
    }

    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    return 0;
}


// ************************************************************************* //
