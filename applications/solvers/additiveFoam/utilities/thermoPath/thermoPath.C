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

#include "thermoPath.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::thermoPath::transportPropertiesName
(
    "transportProperties"
);

const Foam::word Foam::thermoPath::thermoPathName
(
    "thermoPath"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::thermoPath::read
(
    const IOdictionary& transportProperties,
    const Time& runTime
)
{
    if (transportProperties.found(thermoPathName))
    {
        const List<Pair<scalar>> path
        (
            transportProperties.lookup(thermoPathName)
        );

        set(path, transportProperties.name() + ":" + thermoPathName);
    }
    else
    {
        readLegacyFile(runTime);
    }
}


void Foam::thermoPath::readLegacyFile(const Time& runTime)
{
    const fileName pathName
    (
        runTime.rootPath()
       /runTime.globalCaseName()
       /runTime.constant()
       /thermoPathName
    );

    IFstream is(pathName);

    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot read " << thermoPathName
            << " from " << transportPropertiesName
            << " and cannot open legacy file " << pathName
            << exit(FatalError);
    }

    const List<Pair<scalar>> path(is);

    set(path, pathName);
}


void Foam::thermoPath::set
(
    const List<Pair<scalar>>& path,
    const fileName& sourceName
)
{
    sourceName_ = sourceName;

    if (path.size() < 2)
    {
        FatalErrorInFunction
            << "Thermodynamic path " << sourceName_
            << " must contain at least two (T alpha.solid) points"
            << exit(FatalError);
    }

    T_.setSize(path.size());
    alpha_.setSize(path.size());

    forAll(path, i)
    {
        T_[i] = path[i].first();
        alpha_[i] = path[i].second();
    }

    solidus_ = temperatureAtAlpha(1.0);
    liquidus_ = temperatureAtAlpha(0.0);
}


Foam::scalar Foam::thermoPath::temperatureAtAlpha(const scalar alpha) const
{
    forAll(alpha_, i)
    {
        if (mag(alpha_[i] - alpha) <= small)
        {
            return T_[i];
        }
    }

    for (label i=1; i<alpha_.size(); ++i)
    {
        const scalar alpha0 = alpha_[i - 1];
        const scalar alpha1 = alpha_[i];

        if
        (
            ((alpha - alpha0)*(alpha - alpha1) <= 0)
         && (mag(alpha1 - alpha0) > small)
        )
        {
            return
                T_[i - 1]
              + (alpha - alpha0)/(alpha1 - alpha0)*(T_[i] - T_[i - 1]);
        }
    }

    FatalErrorInFunction
        << "Thermodynamic path " << sourceName_
        << " does not bracket alpha.solid = " << alpha
        << ". Add a path point at or across this solid fraction."
        << exit(FatalError);

    return 0;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::thermoPath::thermoPath(const objectRegistry& db)
:
    T_(),
    alpha_(),
    solidus_(great),
    liquidus_(great),
    sourceName_()
{
    if (db.foundObject<IOdictionary>(transportPropertiesName))
    {
        read(db.lookupObject<IOdictionary>(transportPropertiesName), db.time());
    }
    else
    {
        IOdictionary transportProperties
        (
            IOobject
            (
                transportPropertiesName,
                db.time().constant(),
                db,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        read(transportProperties, db.time());
    }
}


Foam::thermoPath::thermoPath
(
    const IOdictionary& transportProperties,
    const Time& runTime
)
:
    T_(),
    alpha_(),
    solidus_(great),
    liquidus_(great),
    sourceName_()
{
    read(transportProperties, runTime);
}

// ************************************************************************* //
