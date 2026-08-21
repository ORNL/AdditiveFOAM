/*---------------------------------------------------------------------------*\
-------------------------------------------------------------------------------
                Copyright (C) 2026 Oak Ridge National Laboratory
-------------------------------------------------------------------------------
Application
    tabulatedProfileInfo

Description
    Report integral and second-moment metrics for an AdditiveFOAM tabulated
    heat-source profile.
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOmanip.H"
#include "IOstreams.H"
#include "tabulatedProfile.H"

using namespace Foam;

int main(int argc, char* argv[])
{
    argList::addNote
    (
        "Report metrics for an AdditiveFOAM tabulated heat-source profile."
    );
    argList::noParallel();
    argList::validArgs.append("tabulated profile");

    argList args(argc, argv);

    const fileName profileFile(args[1]);
    const tabulatedProfile profile(profileFile);

    Info<< setprecision(16)
        << "profileFile: " << profileFile << nl
        << "nx: " << profile.nx() << nl
        << "ny: " << profile.ny() << nl
        << "x0: " << profile.x0() << nl
        << "x1: " << profile.x1() << nl
        << "y0: " << profile.y0() << nl
        << "y1: " << profile.y1() << nl
        << "dx: " << profile.dx() << nl
        << "dy: " << profile.dy() << nl
        << "integral: " << profile.integral() << nl
        << "centroidX: " << profile.centroidX() << nl
        << "centroidY: " << profile.centroidY() << nl
        << "D4SigmaMajor: " << profile.D4SigmaMajor() << nl
        << "D4SigmaMinor: " << profile.D4SigmaMinor() << nl
        << "D4Sigma: " << profile.D4Sigma() << nl
        << "azimuthRadians: " << profile.azimuth() << endl;

    return 0;
}

// ************************************************************************* //
