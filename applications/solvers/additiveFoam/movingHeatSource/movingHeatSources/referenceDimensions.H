/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2026 Oak Ridge National Laboratory
     \\/     M anipulation  |
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

#ifndef referenceDimensions_H
#define referenceDimensions_H

#include "boundBox.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "profileMetrics.H"

namespace Foam
{

class referenceDimensions
{
    enum D4SigmaComponent
    {
        areaEquivalent,
        major,
        minor
    };

    const fvMesh& mesh_;

    word widthReference_;

    word depthReference_;

    D4SigmaComponent component_;

    scalar isotherm_;

    scalar width_;

    scalar depth_;

public:

    referenceDimensions(const dictionary& dict, const fvMesh& mesh);

    void update
    (
        const profileMetrics& metrics,
        const scalar depth,
        const vector& position,
        const boundBox& bounds
    );

    inline scalar width() const
    {
        return width_;
    }

    inline scalar depth() const
    {
        return depth_;
    }

    inline scalar aspectRatio() const
    {
        return 2.0*depth_/width_;
    }

    inline bool isothermDepth() const
    {
        return depthReference_ == "isotherm";
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
