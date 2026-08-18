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

Class
    Foam::tabulatedProfile

Description
    Uniform two-dimensional profile with exact integral and second moments of
    its bilinear interpolant.

SourceFiles
    tabulatedProfile.C

\*---------------------------------------------------------------------------*/

#ifndef tabulatedProfile_H
#define tabulatedProfile_H

#include "fileName.H"
#include "scalarField.H"

namespace Foam
{

class tabulatedProfile
{
    // Private Data

        label nx_;
        label ny_;

        scalar x0_;
        scalar y0_;
        scalar dx_;
        scalar dy_;
        scalar x1_;
        scalar y1_;
        scalar invDx_;
        scalar invDy_;

        scalarField values_;

        scalar integral_;
        scalar centroidX_;
        scalar centroidY_;
        scalar d4SigmaMajor_;
        scalar d4SigmaMinor_;
        scalar d4SigmaEquivalent_;
        scalar azimuth_;

    // Private Member Functions

        void integrate();


public:

    // Constructors

        tabulatedProfile();

        explicit tabulatedProfile(const fileName& profileFile);

    // Member Functions

        void read(const fileName& profileFile);

        inline scalar interpolate(const scalar x, const scalar y) const
        {
            const scalar xp = (x - x0_)*invDx_;
            const scalar yp = (y - y0_)*invDy_;

            if (xp < 0 || xp > nx_ - 1 || yp < 0 || yp > ny_ - 1)
            {
                return 0;
            }

            const label i = min(label(xp), nx_ - 2);
            const label j = min(label(yp), ny_ - 2);
            const scalar tx = xp - i;
            const scalar ty = yp - j;
            const label index = i + nx_*j;

            return
                (1.0 - tx)*(1.0 - ty)*values_[index]
              + tx*(1.0 - ty)*values_[index + 1]
              + (1.0 - tx)*ty*values_[index + nx_]
              + tx*ty*values_[index + nx_ + 1];
        }

        inline label nx() const
        {
            return nx_;
        }

        inline label ny() const
        {
            return ny_;
        }

        inline scalar x0() const
        {
            return x0_;
        }

        inline scalar y0() const
        {
            return y0_;
        }

        inline scalar x1() const
        {
            return x1_;
        }

        inline scalar y1() const
        {
            return y1_;
        }

        inline scalar dx() const
        {
            return dx_;
        }

        inline scalar dy() const
        {
            return dy_;
        }

        inline scalar integral() const
        {
            return integral_;
        }

        inline scalar centroidX() const
        {
            return centroidX_;
        }

        inline scalar centroidY() const
        {
            return centroidY_;
        }

        inline scalar d4SigmaMajor() const
        {
            return d4SigmaMajor_;
        }

        inline scalar d4SigmaMinor() const
        {
            return d4SigmaMinor_;
        }

        inline scalar d4SigmaEquivalent() const
        {
            return d4SigmaEquivalent_;
        }

        inline scalar azimuth() const
        {
            return azimuth_;
        }
};

} // End namespace Foam

#endif

// ************************************************************************* //
