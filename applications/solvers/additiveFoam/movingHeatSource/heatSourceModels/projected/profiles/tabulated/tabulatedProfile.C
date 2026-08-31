/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                Copyright (C) 2026 Oak Ridge National Laboratory
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "tabulatedProfile.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "IFstream.H"

namespace Foam
{
namespace heatSourceProfiles
{
    defineTypeNameAndDebug(tabulated, 0);
    addToRunTimeSelectionTable
    (
        heatSourceProfile,
        tabulated,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatSourceProfiles::tabulated::tabulated()
:
    nx_(0),
    ny_(0),
    x0_(0),
    y0_(0),
    dx_(0),
    dy_(0),
    x1_(0),
    y1_(0),
    invDx_(0),
    invDy_(0),
    values_(0),
    metrics_(),
    bounds_(point::zero, point::zero)
{}

Foam::heatSourceProfiles::tabulated::tabulated(const fileName& profileFile)
:
    tabulated()
{
    read(profileFile);
}


Foam::heatSourceProfiles::tabulated::tabulated
(
    const dictionary& dict,
    const fvMesh& mesh,
    const scalar
)
:
    tabulated()
{
    const fileName fName(dict.lookup("file"));

    const fileName tableFile
    (
        mesh.time().rootPath()
       /mesh.time().globalCaseName()
       /mesh.time().constant()
       /fName
    );

    read(tableFile);

    Info<< "Tabulated profile: integral=" << metrics_.integral()
        << ", centroid=" << metrics_.centroid() << " m"
        << ", D4Sigma (major minor)=" << metrics_.D4Sigma() << " m"
        << ", azimuth=" << metrics_.azimuth() << " rad" << endl;
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //

void Foam::heatSourceProfiles::tabulated::integrate()
{
    // Raw moments M_pq = integral x^p*y^q*I(x,y) dx dy
    scalar M00 = 0;
    scalar M10 = 0;
    scalar M01 = 0;
    scalar M20 = 0;
    scalar M02 = 0;
    scalar M11 = 0;

    for (label j=0; j<ny_ - 1; ++j)
    {
        for (label i=0; i<nx_ - 1; ++i)
        {
            const label id = i + nx_*j;

            const scalar f00 = values_[id];
            const scalar f10 = values_[id + 1];
            const scalar f01 = values_[id + nx_];
            const scalar f11 = values_[id + nx_ + 1];

            const scalar m00 = 0.25*(f00 + f10 + f01 + f11);
            const scalar m10 = (f00 + 2.0*f10 + f01 + 2.0*f11)/12.0;
            const scalar m01 = (f00 + f10 + 2.0*f01 + 2.0*f11)/12.0;
            const scalar m20 = (f00 + 3.0*f10 + f01 + 3.0*f11)/24.0;
            const scalar m02 = (f00 + f10 + 3.0*f01 + 3.0*f11)/24.0;
            const scalar m11 = (f00 + 2.0*f10 + 2.0*f01 + 4.0*f11)/36.0;

            const scalar x = x0_ + i*dx_;
            const scalar y = y0_ + j*dy_;
            const scalar dA = dx_*dy_;

            M00 += dA*m00;
            M10 += dA*(x*m00 + dx_*m10);
            M01 += dA*(y*m00 + dy_*m01);
            M20 +=
                dA*(sqr(x)*m00 + 2.0*x*dx_*m10 + sqr(dx_)*m20);
            M02 +=
                dA*(sqr(y)*m00 + 2.0*y*dy_*m01 + sqr(dy_)*m02);
            M11 +=
                dA
               *(
                    x*y*m00 + x*dy_*m01 + y*dx_*m10
                  + dx_*dy_*m11
                );
        }
    }

    metrics_.reset(M00, M10, M01, M20, M02, M11);
}

// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::heatSourceProfiles::tabulated::read(const fileName& profileFile)
{
    IFstream is(profileFile);

    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot open tabulated profile " << profileFile
            << exit(FatalError);
    }

    is >> nx_ >> ny_;
    is >> x0_ >> y0_;
    is >> dx_ >> dy_;

    if (nx_ < 2 || ny_ < 2 || dx_ <= 0 || dy_ <= 0)
    {
        FatalErrorInFunction
            << "Tabulated profile requires nx, ny >= 2 and positive spacing; "
            << "found nx=" << nx_ << ", ny=" << ny_ << ", dx=" << dx_
            << ", dy=" << dy_ << exit(FatalError);
    }

    x1_ = x0_ + (nx_ - 1)*dx_;
    y1_ = y0_ + (ny_ - 1)*dy_;
    invDx_ = 1.0/dx_;
    invDy_ = 1.0/dy_;

    values_.setSize(nx_*ny_);

    forAll(values_, i)
    {
        is >> values_[i];

        if (!is.good() || values_[i] < 0)
        {
            FatalErrorInFunction
                << "Tabulated profile values must be present and nonnegative; "
                << "invalid value at index " << i << " in " << profileFile
                << exit(FatalError);
        }
    }

    integrate();

    bounds_ =
        boundBox
        (
            point(x0_, y0_, 0),
            point(x1_, y1_, 0)
        );
}

// ************************************************************************* //
