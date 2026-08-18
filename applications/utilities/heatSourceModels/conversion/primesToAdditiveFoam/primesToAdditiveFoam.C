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
#include "tabulatedProfile.H"

using namespace Foam;

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
    const scalar x0 = -0.5 * scalar(nx - 1) * dx;
    const scalar y0 = -0.5 * scalar(ny - 1) * dy;

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
        os.precision(16);

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

    const tabulatedProfile profile(outputFile);

    Info<< "Input: " << inputFile << nl
        << "Output: " << outputFile << nl
        << "Grid: " << nx << " x " << ny << nl
        << "Spacing: " << dx << " x " << dy << " m" << nl
        << "Input integral: " << integral << nl
        << "Normalized integral: " << profile.integral() << nl
        << "D4sigma equivalent: " << profile.d4SigmaEquivalent() << " m" << nl
        << "D4sigma major: " << profile.d4SigmaMajor() << " m" << nl
        << "D4sigma minor: " << profile.d4SigmaMinor() << " m" << nl
        << "Azimuth: " << profile.azimuth() << " rad" << endl;

    return 0;
}

// ************************************************************************* //
