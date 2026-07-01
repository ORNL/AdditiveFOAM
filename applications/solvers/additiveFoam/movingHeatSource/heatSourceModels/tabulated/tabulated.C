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
#include "IFstream.H"

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
    A_(Zero),
    B_(Zero),
    k_(Zero),
    x0_(Zero),
    y0_(Zero),
    dx_(Zero),
    dy_(Zero),
    invDx_(Zero),
    invDy_(Zero),
    nx_(Zero),
    ny_(Zero),
    f_(0),
    Axy_(Zero)
{
    A_ = heatSourceModelCoeffs_.lookup<scalar>("A");
    B_ = heatSourceModelCoeffs_.lookup<scalar>("B");

    const fileName fName(heatSourceModelCoeffs_.lookup("file"));

    const fileName tableFile
    (
        mesh_.time().rootPath()
       /mesh_.time().globalCaseName()
       /mesh_.time().constant()
       /fName
    );

    readTable(tableFile);
    integrateTable();

    const scalar x =
        max
        (
            dimensions_.z()/min(staticDimensions_.x(), staticDimensions_.y()),
            0.001
        );

    const scalar n = min(max(A_*std::log2(x) + B_, 0.0), 9.0);

    k_ = std::pow(2.0, n);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::heatSourceModels::tabulated::readTable(const fileName& tableFile)
{
    IFstream is(tableFile);

    if (!is.good())
    {
        FatalIOErrorInFunction(heatSourceModelCoeffs_)
            << "Cannot find tabulated heat source file " << tableFile
            << exit(FatalIOError);
    }

    is >> nx_ >> ny_;
    is >> x0_ >> y0_;
    is >> dx_ >> dy_;

    invDx_ = 1.0/dx_;
    invDy_ = 1.0/dy_;

    f_.setSize(nx_*ny_);

    forAll(f_, i)
    {
        is >> f_[i];
    }
}


void Foam::heatSourceModels::tabulated::integrateTable()
{
    Axy_ = Zero;

    for (label j=0; j<ny_ - 1; ++j)
    {
        for (label i=0; i<nx_ - 1; ++i)
        {
            const label id = i + nx_*j;

            Axy_ +=
                0.25*dx_*dy_
              * (
                    f_[id]
                  + f_[id + 1]
                  + f_[id + nx_]
                  + f_[id + nx_ + 1]
                );
        }
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline Foam::scalar
Foam::heatSourceModels::tabulated::weight(const vector& d)
{
    const scalar xp = (d.x() - x0_)*invDx_;
    const scalar yp = (d.y() - y0_)*invDy_;

    if (xp < 0.0 || xp > nx_ - 1 || yp < 0.0 || yp > ny_ - 1)
    {
        return 0.0;
    }

    label i = label(xp);
    label j = label(yp);

    i = min(i, nx_ - 2);
    j = min(j, ny_ - 2);

    const scalar tx = xp - i;
    const scalar ty = yp - j;

    const label id = i + nx_*j;

    const scalar f =
        (1.0 - tx)*(1.0 - ty)*f_[id]
      + tx*(1.0 - ty)*f_[id + 1]
      + (1.0 - tx)*ty*f_[id + nx_]
      + tx*ty*f_[id + nx_ + 1];

    const scalar s =
        std::exp(-3.0*std::pow(mag(d.z()/dimensions_.z()), k_));

    return f*s;
}


inline Foam::dimensionedScalar
Foam::heatSourceModels::tabulated::V0()
{
    const scalar x =
        max
        (
            dimensions_.z()/min(staticDimensions_.x(), staticDimensions_.y()),
            0.001
        );

    const scalar n = min(max(A_*std::log2(x) + B_, 0.0), 9.0);

    k_ = std::pow(2.0, n);

    const dimensionedScalar V0
    (
        "V0",
        dimVolume,
        Axy_*dimensions_.z()*Foam::tgamma(1.0/k_)
      / (k_*std::pow(3.0, 1.0/k_))
    );

    return V0;
}


bool Foam::heatSourceModels::tabulated::read()
{
    if (heatSourceModel::read())
    {
        heatSourceModelCoeffs_ = optionalSubDict(type() + "Coeffs");

        heatSourceModelCoeffs_.lookup("A") >> A_;
        heatSourceModelCoeffs_.lookup("B") >> B_;

        const fileName fName(heatSourceModelCoeffs_.lookup("file"));

        const fileName tableFile
        (
            mesh_.time().rootPath()
           /mesh_.time().globalCaseName()
           /mesh_.time().constant()
           /fName
        );

        readTable(tableFile);
        integrateTable();

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
