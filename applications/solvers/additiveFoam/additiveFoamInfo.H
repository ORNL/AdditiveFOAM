/*---------------------------------------------------------------------------*\
-------------------------------------------------------------------------------
                Copyright (C) 2023 Oak Ridge National Laboratory                
-------------------------------------------------------------------------------

Class
    Foam::AdditiveFoamInfo

Description
    Outputs version and build information for AdditiveFOAM.

\*---------------------------------------------------------------------------*/

#ifndef AdditiveFoamInfo_H
#define AdditiveFoamInfo_H

#include "OSspecific.H"

#ifdef GIT_MODULE_ENABLED
    #include "Make/gitInfo.H"
#endif

// Static version
#define ADDITIVEFOAM_VERSION "1.1.0"

namespace Foam
{

class AdditiveFoamInfo
{
public:

    //- Disallow default construction, copy, assignment, and destruction
    AdditiveFoamInfo() = delete;
    ~AdditiveFoamInfo() = delete;
    AdditiveFoamInfo(const AdditiveFoamInfo&) = delete;
    AdditiveFoamInfo& operator=(const AdditiveFoamInfo&) = delete;

    //- Print version/build information to Info stream
    static void write()
    {
        Info<< "AdditiveFOAM Information:" << nl;
        Info<< "Version:   " << ADDITIVEFOAM_VERSION << nl;

        //- Print git describe and SHA1 (optional)
        #ifdef ADDITIVEFOAM_GIT_DESCRIBE
            Info<< "Build:     " << ADDITIVEFOAM_GIT_DESCRIBE << nl;
        #endif

        #ifdef ADDITIVEFOAM_GIT_SHA1
            Info<< "Git SHA1:  " << ADDITIVEFOAM_GIT_SHA1 << nl;
        #endif

        Info<< endl;
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //

