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
    const heatSourceProfiles::tabulated profile(profileFile);
    const profileMetrics& metrics = profile.metrics();

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
        << "integral: " << metrics.integral() << nl
        << "centroidX: " << metrics.centroid().x() << nl
        << "centroidY: " << metrics.centroid().y() << nl
        << "D4Sigma: " << metrics.D4Sigma() << nl
        << "azimuthRadians: " << metrics.azimuth() << endl;

    return 0;
}

// ************************************************************************* //
