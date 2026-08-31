/*---------------------------------------------------------------------------*\
-------------------------------------------------------------------------------
                Copyright (C) 2026 Oak Ridge National Laboratory
-------------------------------------------------------------------------------
Application
    primesToAdditiveFoam

Description
    Convert a PRIMES LaserDiagnosticsSoftware CSV beam-profile export to the
    AdditiveFOAM tabulated heat-source profile format.
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOstreams.H"
#include "IFstream.H"
#include "ISstream.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "scalarField.H"
#include "stringList.H"
#include "HashTable.H"
#include "mathematicalConstants.H"
#include "tabulatedProfile.H"
#include "vector.H"

using namespace Foam;
using Foam::constant::mathematical::pi;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace
{

string trim(const string& s)
{
    size_t first = s.find_first_not_of(" \t\r\n");
    if (string::npos == first) return "";
    size_t last = s.find_last_not_of(" \t\r\n");
    return s.substr(first, (last - first + 1));
}

void stripBOM(string& s)
{
    if
    (
        s.size() >= 3
     && static_cast<unsigned char>(s[0]) == 0xef
     && static_cast<unsigned char>(s[1]) == 0xbb
     && static_cast<unsigned char>(s[2]) == 0xbf
    )
    {
        s.erase(0, 3);
    }
}

stringList split(const string& s, char delimiter)
{
    stringList result;
    size_t start = 0;
    size_t end = s.find(delimiter);

    while (end != string::npos)
    {
        result.append(trim(s.substr(start, end - start)));
        start = end + 1;
        end = s.find(delimiter, start);
    }
    result.append(trim(s.substr(start)));
    return result;
}

string lowercase(string s)
{
    for (char &c : s) c = static_cast<char>(tolower(c));
    return s;
}

scalar parseScalarValue(string s)
{
    s = trim(s);
    size_t commaPos = s.find(',');
    if (commaPos != string::npos) s[commaPos] = '.';

    if (s.empty()) return 0.0;

    IStringStream is(s);
    scalar val;
    is >> val;
    return val;
}

void readPrimesFile
(
    const fileName& inputFile,
    HashTable<string>& metadata,
    List<scalarField>& table,
    label& nx,
    label& ny,
    char& separator
)
{
    IFstream is(inputFile);
    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot open file " << inputFile
            << exit(FatalError);
    }

    separator = ';';
    bool readingPixels = false;
    label rowi = 0;
    string line;

    while (is.getLine(line))
    {
        stripBOM(line);
        line = trim(line);
        if (line.empty()) continue;

        if (line.size() >= 4 && line.substr(0, 4) == "sep=")
        {
            separator = line[4];
            continue;
        }

        stringList fields = split(line, separator);

        if (!readingPixels)
        {
            if (fields.size() > 0 && lowercase(fields[0]) == "pixel")
            {
                if
                (
                    !metadata.found("# Pixel in x")
                 || !metadata.found("# Pixel in y")
                )
                {
                    FatalErrorInFunction
                        << "Grid size metadata not found before 'Pixel' marker"
                        << exit(FatalError);
                }

                nx = label(parseScalarValue(metadata["# Pixel in x"]));
                ny = label(parseScalarValue(metadata["# Pixel in y"]));

                table.setSize(ny);
                forAll(table, i)
                {
                    table[i].setSize(nx, 0.0);
                }

                readingPixels = true;
                continue;
            }

            if (fields.size() >= 3)
            {
                metadata.insert(fields[0], fields[2]);
            }
            else if (fields.size() >= 2)
            {
                metadata.insert(fields[0], fields[1]);
            }
        }
        else
        {
            if (rowi >= ny) break;

            label coli = 0;
            forAll(fields, i)
            {
                if (fields[i].empty()) continue;
                if (coli < nx)
                {
                    table[rowi][coli] = parseScalarValue(fields[i]);
                    coli++;
                }
            }
            rowi++;
        }
    }
}

scalar calculateIntegral(const List<scalarField>& table, scalar dx, scalar dy)
{
    scalar sum = 0;
    for (label j=0; j < table.size()-1; ++j)
    {
        for (label i=0; i < table[j].size()-1; ++i)
        {
            sum += 0.25 *
            (   table[j][i]
              + table[j][i+1]
              + table[j+1][i]
              + table[j+1][i+1]
            );
        }
    }
    return sum * dx * dy;
}

// Calculate the centroid and x/y second-moment radii of the camera pixels
void calculatePixelMoments
(
    const List<scalarField>& table,
    const scalar x0,
    const scalar y0,
    const scalar dx,
    const scalar dy,
    vector& centroid,
    vector& radii
)
{
    scalar M00 = 0;
    scalar M10 = 0;
    scalar M01 = 0;
    scalar M20 = 0;
    scalar M02 = 0;

    forAll(table, j)
    {
        const scalar y = y0 + j*dy;

        forAll(table[j], i)
        {
            const scalar x = x0 + i*dx;
            const scalar f = table[j][i];

            M00 += f;
            M10 += x*f;
            M01 += y*f;
            M20 += sqr(x)*f;
            M02 += sqr(y)*f;
        }
    }

    centroid = vector(M10/M00, M01/M00, 0);
    radii =
        2.0
       *vector
        (
            Foam::sqrt(M20/M00 - sqr(centroid.x())),
            Foam::sqrt(M02/M00 - sqr(centroid.y())),
            0
        );
}

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Convert PRIMES CSV beam profiles to AdditiveFOAM format."
    );
    argList::noParallel();

    argList::validArgs.append("input PRIMES CSV file");
    argList::validArgs.append("output tabulated profile");

    argList args(argc, argv);

    const fileName inputFile = args[1];
    const fileName outputFile = args[2];

    HashTable<string> metadata;
    List<scalarField> table;
    label nx = 0;
    label ny = 0;
    char sep = ';';

    readPrimesFile(inputFile, metadata, table, nx, ny, sep);

    const scalar dx = parseScalarValue(metadata["Pixel pitch x"]) * 1e-6;
    const scalar dy = parseScalarValue(metadata["Pixel pitch y"]) * 1e-6;
    scalar x0 = -0.5 * scalar(nx - 1) * dx;
    scalar y0 = -0.5 * scalar(ny - 1) * dy;

    if (metadata.found("Nullvalue"))
    {
        scalar nullVal = parseScalarValue(metadata["Nullvalue"]);
        forAll(table, i)
        {
            forAll(table[i], j)
            {
                table[i][j] = max(0.0, table[i][j] - nullVal);
            }
        }
    }

    // ROI fill factor = beam diameter/ROI width
    const scalar roiFillFactor = 0.5;
    const List<scalarField> unmaskedTable(table);
    vector centroid;
    vector radii;

    calculatePixelMoments(table, x0, y0, dx, dy, centroid, radii);

    for (label iteration=0; iteration<max(nx, ny); ++iteration)
    {
        const scalar xMin = centroid.x() - radii.x()/roiFillFactor;
        const scalar xMax = centroid.x() + radii.x()/roiFillFactor;
        const scalar yMin = centroid.y() - radii.y()/roiFillFactor;
        const scalar yMax = centroid.y() + radii.y()/roiFillFactor;

        bool changed = false;

        forAll(table, j)
        {
            const scalar y = y0 + j*dy;

            forAll(table[j], i)
            {
                const scalar x = x0 + i*dx;
                const scalar f =
                    x >= xMin && x <= xMax && y >= yMin && y <= yMax
                  ? unmaskedTable[j][i]
                  : 0;

                changed = changed || f != table[j][i];
                table[j][i] = f;
            }
        }

        if (!changed)
        {
            break;
        }

        calculatePixelMoments(table, x0, y0, dx, dy, centroid, radii);
    }

    label minColumn = nx - 1;
    label maxColumn = 0;
    label minRow = ny - 1;
    label maxRow = 0;

    forAll(table, j)
    {
        forAll(table[j], i)
        {
            if (table[j][i] > 0)
            {
                minColumn = min(minColumn, i);
                maxColumn = max(maxColumn, i);
                minRow = min(minRow, j);
                maxRow = max(maxRow, j);
            }
        }
    }

    // Retain one zero-valued point around the profile so cropping does not
    // change its bilinear interpolation or moments.
    minColumn = max(minColumn - 1, 0);
    maxColumn = min(maxColumn + 1, nx - 1);
    minRow = max(minRow - 1, 0);
    maxRow = min(maxRow + 1, ny - 1);

    const label croppedNx = maxColumn - minColumn + 1;
    const label croppedNy = maxRow - minRow + 1;
    List<scalarField> croppedTable(croppedNy);

    forAll(croppedTable, j)
    {
        croppedTable[j].setSize(croppedNx);

        forAll(croppedTable[j], i)
        {
            croppedTable[j][i] = table[minRow + j][minColumn + i];
        }
    }

    table.transfer(croppedTable);
    nx = croppedNx;
    ny = croppedNy;
    x0 += minColumn*dx;
    y0 += minRow*dy;

    scalar integral = calculateIntegral(table, dx, dy);

    if (integral > VSMALL)
    {
        forAll(table, i)
        {
            table[i] /= integral;
        }
    }

    mkDir(outputFile.path());
    {
        OFstream os(outputFile);
        os.precision(12);

        os << nx << " " << ny << nl
           << x0 << " " << y0 << nl
           << dx << " " << dy << nl;

        forAll(table, rowi)
        {
            forAll(table[rowi], coli)
            {
                os << table[rowi][coli] << (coli == nx-1 ? "" : " ");
            }
            os << nl;
        }
    }

    const heatSourceProfiles::tabulated profile(outputFile);
    const profileMetrics& metrics = profile.metrics();
    const scalar micronsToMetres = 1e-6;

    Info<< "Input: " << inputFile << nl
        << "Output: " << outputFile << nl
        << "Grid: " << nx << " x " << ny << nl
        << "Spacing: " << dx << " x " << dy << " m" << nl
        << "ROI fill factor: " << roiFillFactor << nl
        << "ROI integral: " << integral << nl
        << "Normalized integral: " << metrics.integral() << nl
        << "Centroid: " << metrics.centroid() << " m" << nl
        << "D4Sigma (major minor): " << metrics.D4Sigma() << " m" << nl
        << "Azimuth: " << metrics.azimuth() << " rad" << endl;

    if
    (
        metadata.found("Radius a")
     && metadata.found("Radius b")
     && metadata.found("Azimuth angle φ")
    )
    {
        const scalar degreesToRadians = pi/180.0;

        const scalar primesRadiusA =
            parseScalarValue(metadata["Radius a"])*micronsToMetres;

        const scalar primesRadiusB =
            parseScalarValue(metadata["Radius b"])*micronsToMetres;

        const scalar primesPhi =
            parseScalarValue(metadata["Azimuth angle φ"])*degreesToRadians;

        const scalar primesRadiusX =
            Foam::sqrt
            (
                sqr(primesRadiusA*std::cos(primesPhi))
              + sqr(primesRadiusB*std::sin(primesPhi))
            );

        const scalar primesRadiusY =
            Foam::sqrt
            (
                sqr(primesRadiusA*std::sin(primesPhi))
              + sqr(primesRadiusB*std::cos(primesPhi))
            );

        const scalar primesD4Sigma =
            2.0*Foam::sqrt(primesRadiusA*primesRadiusB);

        const scalar majorRadius = 0.5*metrics.D4Sigma().x();
        const scalar minorRadius = 0.5*metrics.D4Sigma().y();
        const scalar theta = metrics.azimuth();

        const scalar radiusX =
            Foam::sqrt
            (
                sqr(majorRadius*std::cos(theta))
              + sqr(minorRadius*std::sin(theta))
            );

        const scalar radiusY =
            Foam::sqrt
            (
                sqr(majorRadius*std::sin(theta))
              + sqr(minorRadius*std::cos(theta))
            );

        scalar radiusA = majorRadius;
        scalar radiusB = minorRadius;
        scalar phi = theta;

        // PRIMES labels the principal axis closest to x as the a-axis.
        if (mag(phi) > 0.25*pi)
        {
            radiusA = minorRadius;
            radiusB = majorRadius;
            phi +=
                phi > 0
              ? -0.5*pi
                : 0.5*pi;
        }

        const scalar D4Sigma =
            Foam::sqrt
            (
                metrics.D4Sigma().x()*metrics.D4Sigma().y()
            );

        Info<< nl
            << "Beam statistics: PRIMES / converted / difference" << nl
            << "Radius a: " << primesRadiusA/micronsToMetres << " / "
            << radiusA/micronsToMetres << " / "
            << (radiusA - primesRadiusA)/micronsToMetres << " um" << nl
            << "Radius b: " << primesRadiusB/micronsToMetres << " / "
            << radiusB/micronsToMetres << " / "
            << (radiusB - primesRadiusB)/micronsToMetres << " um" << nl
            << "Radius x: " << primesRadiusX/micronsToMetres << " / "
            << radiusX/micronsToMetres << " / "
            << (radiusX - primesRadiusX)/micronsToMetres << " um" << nl
            << "Radius y: " << primesRadiusY/micronsToMetres << " / "
            << radiusY/micronsToMetres << " / "
            << (radiusY - primesRadiusY)/micronsToMetres << " um" << nl
            << "D4Sigma: " << primesD4Sigma/micronsToMetres << " / "
            << D4Sigma/micronsToMetres << " / "
            << (D4Sigma - primesD4Sigma)/micronsToMetres
            << " um" << nl
            << "Azimuth: " << primesPhi/degreesToRadians << " / "
            << phi/degreesToRadians << " / "
            << (phi - primesPhi)/degreesToRadians << " deg" << endl;
    }

    if
    (
        metadata.found("Window center position x")
     && metadata.found("Window center position y")
    )
    {
        const scalar windowCenterX =
            parseScalarValue(metadata["Window center position x"]);

        const scalar windowCenterY =
            parseScalarValue(metadata["Window center position y"]);

        Info<< "Centroid in PRIMES coordinates: ("
            << windowCenterX + metrics.centroid().x()/micronsToMetres << ' '
            << windowCenterY + metrics.centroid().y()/micronsToMetres
            << ") um" << endl;
    }

    return 0;
}

// ************************************************************************* //
