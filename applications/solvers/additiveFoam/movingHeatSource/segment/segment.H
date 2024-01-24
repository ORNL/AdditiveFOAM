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
    segment

Description
    Class for properties of a moving heat source
    
    mode:      0 or 1 (0 = line source, 1 = point source)
    position:  position of the heat source center
    power:     power of the heat source
    parameter: (mode = 0: time interval, mode = 1: scan velocity)

SourceFiles
    segment.C

\*---------------------------------------------------------------------------*/

#ifndef segment_H
#define segment_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "point.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class segment Declaration
\*---------------------------------------------------------------------------*/

class segment
{
private:

    // Private data

        scalar mode_;

        point position_;

        scalar power_;

        scalar parameter_;
        
        scalar time_;

public:

    // Constructors

        //- Default construction
        segment();

        //- Construct from space-delimited string
        segment(std::string);


    //- Destructor
    virtual ~segment()
    {}

    // Member Functions

        //- Set time to provided value
        void setTime(scalar time);

        //- Set position to provided value
        void setPosition(vector position);

        scalar mode() const
        {
            return mode_;
        }
        
        point position() const
        {
            return position_;    
        }
                
        scalar power() const
        {
            return power_;
        }
        
        scalar parameter() const
        {
            return parameter_;
        }
        
        scalar time() const
        {
            return time_;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
