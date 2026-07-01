/*---------------------------------------------------------------------------*\
-------------------------------------------------------------------------------
                Copyright (C) 2023-2026 Oak Ridge National Laboratory
-------------------------------------------------------------------------------

Class
    Foam::AdditiveFoamInfo

Description
    Outputs version and build information for AdditiveFOAM.

\*---------------------------------------------------------------------------*/

#ifndef additiveFoamInfo_H
#define additiveFoamInfo_H

#include "messageStream.H"

#ifndef ADDITIVEFOAM_VERSION
    #define ADDITIVEFOAM_VERSION "undefined"
#endif

#ifndef ADDITIVEFOAM_GIT_DESCRIBE
    #define ADDITIVEFOAM_GIT_DESCRIBE "undefined"
#endif

#ifndef ADDITIVEFOAM_GIT_SHA1
    #define ADDITIVEFOAM_GIT_SHA1 "undefined"
#endif

namespace Foam
{

class AdditiveFoamInfo
{
public:

    AdditiveFoamInfo() = delete;

    ~AdditiveFoamInfo() = delete;

    AdditiveFoamInfo(const AdditiveFoamInfo&) = delete;

    AdditiveFoamInfo& operator=(const AdditiveFoamInfo&) = delete;

    static void write()
    {
        Info<< "AdditiveFOAM Information:" << nl
            << "Version:   " << ADDITIVEFOAM_VERSION << nl
            << "Build:     " << ADDITIVEFOAM_GIT_DESCRIBE << nl
            << "Git SHA1:  " << ADDITIVEFOAM_GIT_SHA1 << nl
            << endl;
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
