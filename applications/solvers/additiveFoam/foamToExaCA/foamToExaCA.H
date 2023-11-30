/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                Copyright (C) 2023 Oak Ridge National Laboratory                
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

Class
    Foam::foamToExaca

Description
    A utility for interpolating temperature data to a uniform cartestian grid
    and writes the data in the ExaCA data format.

SourceFiles
    foamToExaca.C

\*---------------------------------------------------------------------------*/

#ifndef foamToExaCA_H
#define foamToExaCA_H
    
#include "IOdictionary.H"
#include "volFields.H"
#include "volPointInterpolation.H"
#include "meshSearch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class foamToExaca Declaration
\*---------------------------------------------------------------------------*/

class foamToExaCA
:
    public IOdictionary
{
private:
    // Private data

        const fvMesh& mesh_;

        const Time& runTime_;

        const volScalarField& T_;

        volPointInterpolation vpi_;

        mutable pointScalarField Tp_;


        List<DynamicList<label>> pointsInCell;

        DynamicList<point> positions;
        DynamicList<scalarField> weights;

        DynamicList<label> compactCells;

        DynamicList<List<scalar>> events;
        DynamicList<List<scalar>> rates;

        Switch    execute_;

        boundBox  box_;

        scalar    iso_;

        scalar    dx_;

public:
    TypeName("foamToExaCA");

    // Construct from temperature field
    foamToExaCA(const volScalarField& T);

    //- Destructor
    virtual ~foamToExaCA()
    {}

    // Member Functions
    virtual void initialize();

    virtual void update(const volScalarField& R);

    virtual void write();

    virtual void mapPoints(const meshSearch& searchEngine);

    virtual void interpolatePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "foamToExaCA/foamToExaCA.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

