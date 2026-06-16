#ifndef readThermoPath_H
#define readThermoPath_H

#include "Time.H"
#include "fvMesh.H"
#include "IFstream.H"
#include "IOdictionary.H"
#include "graph.H"

namespace Foam
{

inline autoPtr<graph> readThermoPath(const fvMesh& mesh)
{
    const Time& runTime = mesh.time();

    const fileName thermoPathFile
    (
        runTime.rootPath()/runTime.globalCaseName()/runTime.constant()
       /"thermoPath"
    );

    IFstream thermoHeaderStream(thermoPathFile);

    if (!thermoHeaderStream.good())
    {
        FatalErrorInFunction
            << "Cannot read thermoPath file " << thermoPathFile
            << exit(FatalError);
    }

    token firstToken(thermoHeaderStream);

    // Check if thermoPath is legacy AdditiveFOAM v1* version (raw list)
    if (firstToken == token::BEGIN_LIST)
    {
        IFstream thermoFile(thermoPathFile);

        return autoPtr<graph>
        (
            new graph("thermo", "T", "alpha1", thermoFile)
        );
    }
    // Otherwise, should be a foamFile
    else if
    (
        firstToken.isWord()
     && firstToken.wordToken() == IOobject::foamFile
    )
    {
        IOdictionary thermoPathDict
        (
            IOobject
            (
                "thermoPath",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        const dictionary& thermoPathSubDict =
            thermoPathDict.subDict("thermoPath");

        if (!thermoPathSubDict.found("values"))
        {
            FatalErrorInFunction
                << "FoamFile-style thermoPath must contain a 'values' entry "
                << "in sub-dictionary thermoPath"
                << exit(FatalError);
        }

        ITstream& valuesStream = thermoPathSubDict.lookup("values");

        return autoPtr<graph>
        (
            new graph("thermo", "T", "alpha1", valuesStream)
        );
    }

    FatalErrorInFunction
        << "Unsupported thermoPath format in " << thermoPathFile
        << ". Expected a legacy list or a FoamFile dictionary."
        << exit(FatalError);

    return autoPtr<graph>();
}

}

#endif

// ************************************************************************* //
