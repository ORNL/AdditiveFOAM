/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

SourceFiles
    createScanPath.C

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "OFstream.H"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Struct Point Definition
\*---------------------------------------------------------------------------*/

struct Point
{
    scalar x;
    scalar y;

    // Default constructor
    Point()
        : x( 0.0 )
        , y( 0.0 )
    {
    }

    // Construct from x and y
    Point( scalar xCoord, scalar yCoord )
        : x( xCoord )
        , y( yCoord )
    {
    }

    // Function to rotate a point around a specified origin by a given angle
    Point rotate( const Point& origin, scalar degrees ) const
    {
        scalar angle = degrees * ( M_PI / 180.0 );

        scalar s = sin( angle );
        scalar c = cos( angle );
        scalar translatedX = x - origin.x;
        scalar translatedY = y - origin.y;
        scalar newX = translatedX * c - translatedY * s + origin.x;
        scalar newY = translatedX * s + translatedY * c + origin.y;
        return Point( newX, newY );
    }

    // Define a custom operator<< to print a Point
    friend std::ostream& operator<<( std::ostream& os, const Point& point )
    {
        os << "(" << point.x << ", " << point.y << ")";
        return os;
    }
};

// Function to calculate the distance between two points
scalar distance( const Point& p1, const Point& p2 )
{
    scalar dx = p1.x - p2.x;
    scalar dy = p1.y - p2.y;
    return std::sqrt( dx * dx + dy * dy );
}


/*---------------------------------------------------------------------------*\
                            Struct Line Definition
\*---------------------------------------------------------------------------*/

struct Line
{
    Point start;
    Point end;

    // Default constructor
    Line()
        : start( Point() )
        , end( Point() )
    {
    }

    // Construct from two points
    Line( Point startPoint, Point endPoint )
        : start( startPoint )
        , end( endPoint )
    {
    }

    // Function to rotate the line around a specified origin by a given angle
    void rotate( const Point& origin, scalar angle )
    {
        start = start.rotate( origin, angle );
        end = end.rotate( origin, angle );
    }

    bool isFinite()
    {
        return ( !std::isnan( start.x ) && !std::isnan( start.y ) &&
                 !std::isnan( end.x ) && !std::isnan( end.y ) );
    }

    // Define a custom operator<< to print a Line
    friend std::ostream& operator<<( std::ostream& os, const Line& line )
    {
        os << "(" << line.start.x << ", " << line.start.y << "), "
           << "(" << line.end.x << ", " << line.end.y << ")";
        return os;
    }
};


/*---------------------------------------------------------------------------*\
                            Struct BoundBox Definition
\*---------------------------------------------------------------------------*/

struct BoundBox
{
    Point minPoint;
    Point maxPoint;
    Point midPoint;

    std::vector<Line> edges;

    // Construct from two points
    BoundBox( Point minP, Point maxP )
        : minPoint( minP )
        , maxPoint( maxP )
        , midPoint( ( minP.x + maxP.x ) / 2.0, ( minP.y + maxP.y ) / 2.0 )
    {
        // left, right, top, bottom
        edges = { { Point( minPoint.x, minPoint.y ),
                    Point( minPoint.x, maxPoint.y ) },
                  { Point( maxPoint.x, minPoint.y ),
                    Point( maxPoint.x, maxPoint.y ) },
                  { Point( minPoint.x, maxPoint.y ),
                    Point( maxPoint.x, maxPoint.y ) },
                  { Point( minPoint.x, minPoint.y ),
                    Point( maxPoint.x, minPoint.y ) } };
    }

    // Function to check if a point is inside the bounding box
    bool isInside( Point p )
    {
        return p.x >= minPoint.x && p.x <= maxPoint.x && p.y >= minPoint.y &&
               p.y <= maxPoint.y;
    }

    Line cropLine( const Line& line )
    {
        std::vector<Point> intersections;

        for ( const Line& edge : edges )
        {
            Point intersection = intersect( edge, line );

            if ( !std::isnan( intersection.x ) &&
                 !std::isnan( intersection.y ) )
            {
                intersections.push_back( intersection );
            }
        }

        if ( intersections.size() == 0 )
        {
            return Line( Point( NAN, NAN ), Point( NAN, NAN ) );
        }

        // Sort intersection points based on position along the original line
        std::sort( intersections.begin(), intersections.end(),
                   [&line]( const Point& p1, const Point& p2 )
                   {
                       scalar d1 = distance( line.start, p1 );
                       scalar d2 = distance( line.start, p2 );
                       return d1 < d2;
                   } );

        Line intersectedLine( intersections.front(), intersections.back() );

        return intersectedLine;
    }

    // Function to find the intersection point between two lines
    Point intersect( const Line& line1, const Line& line2 )
    {
        Point intersection( std::numeric_limits<scalar>::quiet_NaN(),
                            std::numeric_limits<scalar>::quiet_NaN() );

        // components of first line
        scalar x1 = line1.start.x;
        scalar y1 = line1.start.y;
        scalar x2 = line1.end.x;
        scalar y2 = line1.end.y;

        // components of second line
        scalar x3 = line2.start.x;
        scalar y3 = line2.start.y;
        scalar x4 = line2.end.x;
        scalar y4 = line2.end.y;

        scalar denominator =
            ( x1 - x2 ) * ( y3 - y4 ) - ( y1 - y2 ) * ( x3 - x4 );

        // Lines are parallel or colinear, no intersection
        if ( denominator == 0 )
        {
            return intersection;
        }

        scalar t = ( ( x1 - x3 ) * ( y3 - y4 ) - ( y1 - y3 ) * ( x3 - x4 ) ) /
                   denominator;

        scalar u = -( ( x1 - x2 ) * ( y1 - y3 ) - ( y1 - y2 ) * ( x1 - x3 ) ) /
                   denominator;

        if ( t >= 0 && t <= 1 && u >= 0 && u <= 1 )
        {
            intersection.x = x1 + t * ( x2 - x1 );
            intersection.y = y1 + t * ( y2 - y1 );
            return intersection;
        }

        return intersection;
    }
};


/*---------------------------------------------------------------------------*\
                            Struct Path Definition
\*---------------------------------------------------------------------------*/

struct Path
{
    std::vector<Line> lines;
    scalar power;
    scalar speed;
    scalar dwellTime;

    // Construct the path from a bounding box and hatch spacing

    // Default constructor
    Path() {}

    // Constructor from components
    Path(const BoundBox& bbox, const scalar step, const scalar angle)
    {
        createPath(bbox, step, angle);
    }
    
    // Function to create path from components
    void createPath(BoundBox bbox, const scalar step, const scalar angle)
    {
        label n = countLines( bbox, step );

        // create a pad of infinitely long, equally parallel lines
        std::vector<Line> pathLines;

        const scalar great = 1e10;

        // Create lines in the negative direction, excluding the midpoint line
        for ( label i = n - 1; i > 0; --i )
        {
            scalar height = bbox.midPoint.y - i * step;
            Line currentLine( Point( -great, height ), Point( great, height ) );
            pathLines.push_back( currentLine );
        }

        // Create lines in the positive direction, including the midpoint
        for ( label i = 0; i < n; ++i )
        {
            scalar height = bbox.midPoint.y + i * step;
            Line currentLine( Point( -great, height ), Point( great, height ) );
            pathLines.push_back( currentLine );
        }

        // apply rotation and cropping to the initial lines
        for ( const Line& pathLine : pathLines )
        {
            // Rotate the endpoints by the specified angle
            Line rotatedLine = pathLine;
            rotatedLine.rotate( bbox.midPoint, angle );

            Line croppedLine = bbox.cropLine( rotatedLine );

            if ( croppedLine.isFinite() )
            {
                lines.push_back( croppedLine );
            }
        }
    }

    // Function to find the number of scan vectors in the bounding box
    label countLines( const BoundBox& bbox, scalar step ) const
    {
        label nX = 0;
        label nY = 0;

        for ( scalar x = bbox.minPoint.x; x <= bbox.maxPoint.x; x += step )
        {
            nX++;
        }

        for ( scalar y = bbox.minPoint.y; y <= bbox.maxPoint.y; y += step )
        {
            nY++;
        }

        return ( nX > nY ) ? nX : nY;
    }

    void write( const std::string& filename,
                const bool bi_direction = true ) const
    {
        std::ofstream file( filename );

        file << "Mode\tX(m)\tY(m)\tZ(m)\tPower(W)\ttParam" << std::endl;

        for ( size_t i = 0; i < lines.size(); ++i )
        {
            const Line& line = lines[i];

            Point first = line.start;
            Point second = line.end;

            // reverse the odd lines for bi_directional
            if ( bi_direction && i % 2 == 1 )
            {
                first = line.end;
                second = line.start;
            }

            // hatch (with skywrite)
            if ( i == 0 )
            {
                // no initial dwell
                file << "1\t" << first.x << "\t" << first.y << "\t0\t0\t0\n";
            }
            else
            {
                file << "1\t" << first.x << "\t" << first.y << "\t0\t" << 0
                     << "\t" << dwellTime << "\n";
            }

            // raster
            file << "0\t" << second.x << "\t" << second.y << "\t0\t" << power
                 << "\t" << speed << "\n";
        }

        file.close();
    }
};

} // End namespace Foam

// ************************************************************************* //
