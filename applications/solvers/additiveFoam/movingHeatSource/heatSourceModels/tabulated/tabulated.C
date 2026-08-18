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

#include "tabulated.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatSourceModels
{
    defineTypeNameAndDebug(tabulated, 0);
    addToRunTimeSelectionTable(heatSourceModel, tabulated, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatSourceModels::tabulated::tabulated
(
    const word& sourceName,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    heatSourceModel(typeName, sourceName, dict, mesh),
    mesh_(mesh),
    nSlope_(Zero),
    nIntercept_(Zero),
    axialExponent_(1),
    minimumDepth_(Zero),
    profile_()
{
    nSlope_ =
        heatSourceModelCoeffs_.lookup<scalar>("nSlope");

    nIntercept_ =
        heatSourceModelCoeffs_.lookup<scalar>("nIntercept");

    minimumDepth_ =
        heatSourceModelCoeffs_.lookup<scalar>("minimumDepth");

    if (minimumDepth_ <= 0)
    {
        FatalIOErrorInFunction(heatSourceModelCoeffs_)
            << "minimumDepth must be greater than zero"
            << exit(FatalIOError);
    }

    const fileName fName(heatSourceModelCoeffs_.lookup("file"));

    const fileName tableFile
    (
        mesh_.time().rootPath()
       /mesh_.time().globalCaseName()
       /mesh_.time().constant()
       /fName
    );

    profile_.read(tableFile);

    dimensions_ =
        vector
        (
            0.5*profile_.d4SigmaEquivalent(),
            0.5*profile_.d4SigmaEquivalent(),
            minimumDepth_
        );

    staticDimensions_ = dimensions_;
    measuredDepth_ = transient_ ? 0 : minimumDepth_;

    updateAxialState();

    Info<< "Tabulated profile: integral=" << profile_.integral()
        << ", centroid=(" << profile_.centroidX()
        << ' ' << profile_.centroidY() << ") m"
        << ", D4sigma equivalent/major/minor="
        << profile_.d4SigmaEquivalent() << '/'
        << profile_.d4SigmaMajor() << '/'
        << profile_.d4SigmaMinor() << " m"
        << ", azimuth=" << profile_.azimuth() << " rad" << endl;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::heatSourceModels::tabulated::updateAxialState()
{
    const scalar sourceAspectRatio =
        max
        (
            2.0*dimensions_.z()/profile_.d4SigmaEquivalent(),
            0.001
        );

    const scalar log2AxialExponent =
        min
        (
            max
            (
                nSlope_*std::log2(sourceAspectRatio)
              + nIntercept_,
                0.0
            ),
            9.0
        );

    axialExponent_ = std::pow(2.0, log2AxialExponent);

    static const scalar cutoffTolerance = 1.0e-3;

    const scalar cutoffDepth =
        dimensions_.z()
       *std::pow
        (
            -std::log(cutoffTolerance)/3.0,
            1.0/axialExponent_
        );

    sourceLowerBound_ =
        vector(profile_.x0(), profile_.y0(), -cutoffDepth);

    sourceUpperBound_ =
        vector(profile_.x1(), profile_.y1(), 0);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline Foam::scalar
Foam::heatSourceModels::tabulated::weight(const vector& d)
{
    if (d.z() > small)
    {
        return 0;
    }

    const scalar planarWeight = profile_.interpolate(d.x(), d.y());

    const scalar axialWeight =
        std::exp
        (
            -3.0
           *std::pow(max(-d.z(), scalar(0))/dimensions_.z(), axialExponent_)
        );

    return planarWeight*axialWeight;
}


inline Foam::dimensionedScalar
Foam::heatSourceModels::tabulated::V0()
{
    updateAxialState();

    const dimensionedScalar V0
    (
        "V0",
        dimVolume,
        profile_.integral()*dimensions_.z()
       *Foam::tgamma(1.0/axialExponent_)
       /(axialExponent_*std::pow(3.0, 1.0/axialExponent_))
    );

    return V0;
}


bool Foam::heatSourceModels::tabulated::read()
{
    if (heatSourceModel::read())
    {
        heatSourceModelCoeffs_ = optionalSubDict(type() + "Coeffs");

        heatSourceModelCoeffs_.lookup("nSlope") >> nSlope_;
        heatSourceModelCoeffs_.lookup("nIntercept")
            >> nIntercept_;
        heatSourceModelCoeffs_.lookup("minimumDepth") >> minimumDepth_;

        if (minimumDepth_ <= 0)
        {
            FatalIOErrorInFunction(heatSourceModelCoeffs_)
                << "minimumDepth must be greater than zero"
                << exit(FatalIOError);
        }

        const fileName fName(heatSourceModelCoeffs_.lookup("file"));

        const fileName tableFile
        (
            mesh_.time().rootPath()
           /mesh_.time().globalCaseName()
           /mesh_.time().constant()
           /fName
        );

        profile_.read(tableFile);

        const scalar effectiveDepth =
            transient_ ? max(minimumDepth_, measuredDepth_) : minimumDepth_;

        dimensions_ =
            vector
            (
                0.5*profile_.d4SigmaEquivalent(),
                0.5*profile_.d4SigmaEquivalent(),
                effectiveDepth
            );

        staticDimensions_ =
            vector(dimensions_.x(), dimensions_.y(), minimumDepth_);

        if (!transient_)
        {
            measuredDepth_ = minimumDepth_;
        }

        updateAxialState();

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
