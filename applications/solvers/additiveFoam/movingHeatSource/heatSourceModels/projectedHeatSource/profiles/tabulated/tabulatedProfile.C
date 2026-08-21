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
    defineTypeNameAndDebug(tabulatedProfile, 0);
    addToRunTimeSelectionTable
    (
        heatSourceProfile,
        tabulatedProfile,
        dictionary
    );
}

// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tabulatedProfile::tabulatedProfile()
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
    integral_(0),
    centroidX_(0),
    centroidY_(0),
    D4SigmaMajor_(0),
    D4SigmaMinor_(0),
    D4Sigma_(0),
    azimuth_(0),
    profileBb_(point::zero, point::zero)
{}

Foam::tabulatedProfile::tabulatedProfile(const fileName& profileFile)
:
    tabulatedProfile()
{
    read(profileFile);
}


Foam::tabulatedProfile::tabulatedProfile
(
    const dictionary& dict,
    const fvMesh& mesh,
    const scalar
)
:
    tabulatedProfile()
{
    const dictionary& coeffs = dict.subDict(typeName + "Coeffs");

    const fileName fName(coeffs.lookup("file"));

    const fileName tableFile
    (
        mesh.time().rootPath()
       /mesh.time().globalCaseName()
       /mesh.time().constant()
       /fName
    );

    read(tableFile);

    Info<< "Tabulated profile: integral=" << integral_
        << ", centroid=(" << centroidX_ << ' ' << centroidY_ << ") m"
        << ", D4Sigma/major/minor=" << D4Sigma_ << '/'
        << D4SigmaMajor_ << '/' << D4SigmaMinor_ << " m"
        << ", azimuth=" << azimuth_ << " rad" << endl;
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //

void Foam::tabulatedProfile::integrate()
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
            M20 += dA*(sqr(x)*m00 + 2.0*x*dx_*m10 + sqr(dx_)*m20);
            M02 += dA*(sqr(y)*m00 + 2.0*y*dy_*m01 + sqr(dy_)*m02);
            M11 += dA*(x*y*m00 + x*dy_*m01 + y*dx_*m10 + dx_*dy_*m11);
        }
    }

    integral_ = M00;
    centroidX_ = M10/M00;
    centroidY_ = M01/M00;

    const scalar varianceX = M20/M00 - sqr(centroidX_);

    const scalar varianceY = M02/M00 - sqr(centroidY_);

    const scalar covariance = M11/M00 - centroidX_*centroidY_;

    const scalar meanVariance = 0.5*(varianceX + varianceY);
    const scalar delta =
        Foam::sqrt(0.25*sqr(varianceX - varianceY) + sqr(covariance));

    const scalar majorVariance = meanVariance + delta;
    const scalar minorVariance = meanVariance - delta;

    D4SigmaMajor_ = 4.0*Foam::sqrt(majorVariance);
    D4SigmaMinor_ = 4.0*Foam::sqrt(minorVariance);
    D4Sigma_ = Foam::sqrt(D4SigmaMajor_*D4SigmaMinor_);

    azimuth_ =
        delta > small
      ? 0.5*Foam::atan2(2.0*covariance, varianceX - varianceY)
      : 0.0;
}

// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::tabulatedProfile::read(const fileName& profileFile)
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

    x1_ = x0_ + (nx_ - 1)*dx_;
    y1_ = y0_ + (ny_ - 1)*dy_;
    invDx_ = 1.0/dx_;
    invDy_ = 1.0/dy_;

    values_.setSize(nx_*ny_);

    forAll(values_, i)
    {
        is >> values_[i];
    }

    integrate();

    profileBb_ =
        boundBox
        (
            point(x0_, y0_, 0),
            point(x1_, y1_, 0)
        );
}

// ************************************************************************* //
