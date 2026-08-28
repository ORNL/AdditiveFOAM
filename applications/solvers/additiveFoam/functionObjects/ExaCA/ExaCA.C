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

\*---------------------------------------------------------------------------*/

#include "ExaCA.H"
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "labelVector.H"
#include "pointMVCWeight.H"
#include <mpi.h>
#include <array>
#include <adios2.h>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(ExaCA, 0);
    
    addToRunTimeSelectionTable
    (
        functionObject,
        ExaCA,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::ExaCA::ExaCA
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    T_(mesh_.lookupObject<VolField<scalar>>("T")),
    vpi_(mesh_),
    Tp_
    (
        IOobject
        (
            "Tp_",
            runTime.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        vpi_.interpolate(T_)
    ),
    searchEngine_(mesh_, polyMesh::CELL_TETS)
{
    read(dict);
    
    setOverlapCells();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::ExaCA::~ExaCA()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::ExaCA::read(const dictionary& dict)
{
    box_ = dict.lookup("box");
    
    isoValue_ = dict.lookup<scalar>("isoValue");
    
    dx_  = dict.lookup<scalar>("dx");

    return true;
}

void Foam::functionObjects::ExaCA::setOverlapCells()
{
    // create a compact cell-stencil using the overlap sub-space
    const pointField& points = mesh_.points();

    boundBox procBb(points);

    const vector extend = 1e-10*vector::one;

    if (procBb.overlaps(box_))
    {
        forAll(mesh_.cells(), celli)
        {
            boundBox cellBb(point::max, point::min);

            const labelList& vertices = mesh_.cellPoints()[celli];

            forAll(vertices, i)
            {
                cellBb.min() = min(cellBb.min(), points[vertices[i]] - extend);
                cellBb.max() = max(cellBb.max(), points[vertices[i]] + extend);
            }

            if (cellBb.overlaps(box_))
            {
                overlapCells.append(celli);
            }
        }
    }

    overlapCells.shrink();
}


Foam::wordList Foam::functionObjects::ExaCA::fields() const
{
    return wordList::null();
}


bool Foam::functionObjects::ExaCA::execute()
{    
    const pointScalarField Tp0_("Tp0_", Tp_);

    Tp_ = vpi_.interpolate(T_);

    // capture events at interface cells: {cell id, time, vertex temperatures}
    const scalar t_ = mesh_.time().value();
    const scalar t0_ = t_ - mesh_.time().deltaTValue();
        
    forAll(overlapCells, i)
    {
        label celli = overlapCells[i];

        const labelList& vertices = mesh_.cellPoints()[celli];

        label c0 = 0;
        label c1 = 0;

        forAll(vertices, pI)
        {
            if (Tp0_[vertices[pI]] >= isoValue_)
            {
                c0++;
            }

            if (Tp_[vertices[pI]] >= isoValue_)
            {
                c1++;
            }
        }
        
        const label n = vertices.size();

        // overshoot correction: interface jumped this cell during time step
        if ( !(c0 % n) && !(c1 % n) )
        {
            if  (c0 != c1)
            {
                c0 = 0;
                c1 = 1;
            }
        }
        
        // capture solidification events
        c0 %= n;
        c1 %= n;
        
        if (c0 || c1)
        {
            List<scalar> event(n + 2);

            event[0] = celli;
        
            // add previous event
            if (c1 && !c0)
            {
                event[1] = t0_;

                forAll(vertices, pI)
                {
                    event[pI + 2] = Tp0_[vertices[pI]];
                }

                events.append(event);
            }

            // add current event
            event[1] = t_;

            forAll(vertices, pI)
            {
                event[pI + 2] = Tp_[vertices[pI]];
            }

            events.append(event);
        }
    }
    
    return true;
}

void Foam::functionObjects::ExaCA::mapPoints()
{
    if (events.size() == 0)
    {
        return;
    }

    // find event sub-space before constructing interpolant weights
    const pointField& points = mesh_.points();

    const vector extend = 1e-10*vector::one;

    boundBox eventBb(point::max, point::min);

    for (label i = 1; i < events.size(); i++)
    {
        const label celli = events[i][0];

        if (celli == events[i - 1][0])
        {
            continue;
        }

        const labelList& vertices = mesh_.cellPoints()[celli];

        forAll(vertices, i)
        {
            eventBb.min() = min(eventBb.min(), points[vertices[i]] - extend);
            eventBb.max() = max(eventBb.max(), points[vertices[i]] + extend);
        }
    }

    // find interpolant weigthts for each point
    label pI = 0;
    label seedi = events[0][0];

    pointsInCell.setSize(mesh_.nCells());

    const labelVector nPoints(vector::one + box_.span() / dx_);

    Info << "starting point loop" << endl;

    for (label k=0; k < nPoints.z(); ++k)
    {
        for (label j=0; j < nPoints.y(); ++j)
        {
            for (label i=0; i < nPoints.x(); ++i)
            {
                const point pt = box_.max() - vector(i, j, k)*dx_;

                if (eventBb.contains(pt))
                {
                    // shift point during search to avoid edges in pointMVC
                    const point spt = pt - vector::one*1e-10;

                    label celli = searchEngine_.findCell(spt, seedi, true);

                    if (celli != -1)
                    {
                        positions.append(pt);
                  
                        pointMVCWeight cpw(mesh_, spt, celli);

                        weights.append(cpw.weights());

                        pointsInCell[celli].append(pI);

                        pI++;
                    }

                    seedi = celli;
                }
            }
        }
    }

    positions.shrink();

    weights.shrink();

    for (auto& pic : pointsInCell)
    {
        pic.shrink();
    }
}

void Foam::functionObjects::ExaCA::interpolate()
{    
    // Initialize ADIOS
    //adios2::ADIOS adios(MPI_COMM_WORLD);//Pstream::worldComm);

    // Get number of processes
    // Test print
    //const std::string greeting = "Hello World from MPI rank" + std::to_string(Pstream::myProcNo()); //mpi_rank);
    //std::cout << "Writing ADIOS2 data to file(s)" << std::endl;
    //writer(adios, greeting);

    // All MPI ranks must call ADIOS routines even if there are no local events
    // if (events.size() == 0)
    // {
    //     return;
    // }
    
    // format events in ExaCA reduced data format
    DynamicList<List<scalar>> data;
    List<scalar> tm;
    if (events.size() > 0)
        tm.setSize(pointsInCell[events[0][0]].size(), events[0][1]);

    for (label i = 1; i < events.size(); i++)
    {
        const List<scalar> prevEvent = events[i - 1];

        const List<scalar> currEvent = events[i];

        label celli = currEvent[0];

        if (currEvent[0] != prevEvent[0])
        {
            tm.setSize(pointsInCell[celli].size(), currEvent[1]);
            continue;
        }

        scalar prevTime = prevEvent[1];

        scalar currTime = currEvent[1];

        List<scalar> psi0(prevEvent.size() - 2);
        List<scalar> psi(currEvent.size() - 2);

        for (int j = 0; j < psi.size(); j++)
        {
            psi0[j] = prevEvent[j + 2];
            psi[j]  = currEvent[j + 2];
        }

        int p = 0;
        for (const label& pointi : pointsInCell[celli])
        {
            scalar tp0 = Zero;
            scalar tp  = Zero;

            List<scalar> w = weights[pointi];
            
            for (int j = 0; j < psi.size(); j++)
            {
                tp0 += w[j]*psi0[j];
                tp  += w[j]*psi[j];
            }

            if ((tp <= isoValue_) && (tp0 > isoValue_))
            {
                const point& pt = positions[pointi];

                scalar m_ = min(max((isoValue_ - tp0)/(tp - tp0), 0), 1);

                data.append
                (
                    {
                        pt[0],
                        pt[1],
                        pt[2],
                        tm[p],
                        prevTime + m_*(currTime - prevTime),
                        (tp0 - tp) / (currTime - prevTime)
                    }
                );
            }
            else if ((tp > isoValue_) && (tp0 <= isoValue_))
            {
                scalar m_ = min(max((isoValue_ - tp0)/(tp - tp0), 0), 1);

                tm[p] = prevTime + m_*(currTime - prevTime);
            }

            p++;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Number of events to write on the local MPI rank
    unsigned long num_events_local = data.size();

    // Local min/max for x, y, z and number of cells per direction
    std::vector<double> min_bounds(3), max_bounds(3);
    unsigned long nx, ny, nz;
    if (num_events_local > 0)
    {
	std::vector<double> x_values_arr(num_events_local), y_values_arr(num_events_local), z_values_arr(num_events_local);   
        for(unsigned long i=0; i<num_events_local; i++)
 	{
           x_values_arr[i] = data[i][0];
           y_values_arr[i] = data[i][1];
	   z_values_arr[i] = data[i][2];
	}
        min_bounds[0] = *std::min_element(x_values_arr.begin(), x_values_arr.end());
        min_bounds[1] = *std::min_element(y_values_arr.begin(), y_values_arr.end());
        min_bounds[2] = *std::min_element(z_values_arr.begin(), z_values_arr.end());
        max_bounds[0] = *std::max_element(x_values_arr.begin(), x_values_arr.end());
        max_bounds[1] = *std::max_element(y_values_arr.begin(), y_values_arr.end());
        max_bounds[2] = *std::max_element(z_values_arr.begin(), z_values_arr.end());
	nx = std::lround((max_bounds[0] - min_bounds[0]) / dx_) + 1;
    	ny = std::lround((max_bounds[1] - min_bounds[1]) / dx_) + 1;
   	nz = std::lround((max_bounds[2] - min_bounds[2]) / dx_) + 1;
    }
    else
    {
        for (int k=0; k<3; k++)
        {
           min_bounds[k] = std::numeric_limits<double>::max();
           max_bounds[k] = std::numeric_limits<double>::lowest();
        }
	nx = 0;
	ny = 0;
	nz = 0;
    }
    // Number of cells on this local MPI rank
    unsigned long num_cells = nx * ny * nz;
    std::cout << "Num cells " << num_cells << std::endl;
    // Assign each melting/solidification event to a remelt_event_num, and store the total events per cell
    std::vector<unsigned long> remelt_count(num_cells, 0), remelt_event_num(num_cells);
    // Store the x, y, and z location (in cell units) of each event
    std::vector<int> x_cell(num_cells), y_cell(num_cells), z_cell(num_cells);
    for(unsigned long i=0; i<num_events_local; i++)
    {
        double x_position = data[i][0];
	x_cell[i] = std::lround((x_position - min_bounds[0]) / dx_);
	double y_position = data[i][1];
        y_cell[i] = std::lround((y_position - min_bounds[1]) / dx_);
        double z_position = data[i][2];
        z_cell[i] = std::lround((z_position - min_bounds[2]) / dx_);
	unsigned long index1d = z_cell[i] * nx * ny + y_cell[i] * nx + x_cell[i];
	remelt_event_num[index1d] = remelt_count[index1d];
	remelt_count[index1d]++;
    }
    // Max number of times a location undergoes remelting
    unsigned long max_remelt_count;
    if (num_events_local > 0) 
    	max_remelt_count = *std::max_element(remelt_count.begin(), remelt_count.end());
    else
	max_remelt_count = 0;
    // Need the global max for the total array size
    unsigned long max_remelt_count_global;
    MPI_Reduce(&max_remelt_count, &max_remelt_count_global, 1, MPI_UNSIGNED_LONG, MPI_MAX, 0, MPI_COMM_WORLD);

    // Array size for writing
    unsigned long arr_size = num_cells * max_remelt_count_global;
    std::cout << "Arr size " << arr_size << std::endl;
    // tm, tl, cr are the vectors to be writen to the file
    std::vector<double> tm_write(arr_size, 0.0), tl_write(arr_size, 0.0), cr_write(arr_size, 0.0);
    for(unsigned long i=0; i<num_events_local; i++)
    {
	unsigned long index1d = z_cell[i] * nx * ny + y_cell[i] * nx + x_cell[i];    
        unsigned long remelt_event_num_cell = remelt_event_num[index1d];
	// 1D index of the 4D array (max_num_remelt_events * size z * size y * size x)
	unsigned long index1d_full = remelt_event_num_cell * nx * ny * nz + z_cell[i] * nx * ny + y_cell[i] * nx + x_cell[i];
	tm_write[index1d_full] = data[i][3];
        tl_write[index1d_full] = data[i][4];
        cr_write[index1d_full] = data[i][5];
    }

    // Global min/max for x, y, z
    std::vector<double> min_bounds_global(3), max_bounds_global(3);
    MPI_Allreduce(min_bounds.data(), min_bounds_global.data(), 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(max_bounds.data(), max_bounds_global.data(), 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    std::cout << "Global origin: " << min_bounds_global[0] << ", " << min_bounds_global[1] << ", " << min_bounds_global[2] << std::endl;
    // Init adios writer
    adios2::ADIOS adios(MPI_COMM_WORLD);
    adios2::IO io = adios.DeclareIO("WriteArrays");

    // Global domain size
    unsigned long nx_global = std::lround((max_bounds_global[0] - min_bounds_global[0]) / dx_) + 1;
    unsigned long ny_global = std::lround((max_bounds_global[1] - min_bounds_global[1]) / dx_) + 1;
    unsigned long nz_global = std::lround((max_bounds_global[2] - min_bounds_global[2]) / dx_) + 1;
    // Integer offset from global origin
    unsigned long x_offset = std::lround((min_bounds[0] - min_bounds_global[0]) / dx_);
    unsigned long y_offset = std::lround((min_bounds[1] - min_bounds_global[1]) / dx_);
    unsigned long z_offset = std::lround((min_bounds[2] - min_bounds_global[2]) / dx_);

    std::cout << "Min/max bounds for x are " << min_bounds[0] << " and " << max_bounds[0] << " x off " << x_offset << " y off " << y_offset << std::endl;

    std::cout << "X = " << min_bounds[0] << " through " << max_bounds[0] << ", Y = " << min_bounds[1] << " through " << max_bounds[1] << ", Z = " << min_bounds[2] << " through " << max_bounds[2] << " offsets are " << x_offset << "," << y_offset << "," << z_offset << std::endl;

    // Store attributes of data in adios structs - shape is global, start and count are local
    const adios2::Dims shape{max_remelt_count_global, nz_global, ny_global, nx_global};
    const adios2::Dims start{0, z_offset, y_offset, x_offset};
    const adios2::Dims count{max_remelt_count_global, nz, ny, nx};
    // Data attributes are global - defined the same on all ranks
    std::array<size_t, 4> extents = {max_remelt_count_global, nz_global, ny_global, nx_global};
    std::array<double, 4> spacing = {1.0, dx_, dx_, dx_};
    std::array<double, 4> origin  = {0.0, min_bounds_global[2], min_bounds_global[1], min_bounds_global[0]};
    io.DefineAttribute<size_t>("extents", extents.data(), extents.size());
    io.DefineAttribute<double>("spacing", spacing.data(), spacing.size());
    io.DefineAttribute<double>("origin", origin.data(), origin.size());
    io.DefineAttribute<std::string>("axis_order", "N,Z,Y,X");

    // Compression
    //auto zfpOp = io.DefineOperator("ZFPCompressor", adios2::ops::ZFP);

    // Define arrays for writing
    auto tm_var = io.DefineVariable<double>("tm", shape, start, count);
    auto tl_var = io.DefineVariable<double>("tl", shape, start, count);
    auto cr_var = io.DefineVariable<double>("cr", shape, start, count);

    // Compress arrays
    //tm_var.AddOperation(zfpOp, {{"accuracy", "0.00001"}});
    //tl_var.AddOperation(zfpOp, {{"accuracy", "0.00001"}});
    //cr_var.AddOperation(zfpOp, {{"accuracy", "0.00001"}});

    // Write arrays
    adios2::Engine writer = io.Open("adios_output.bp", adios2::Mode::Write);
    writer.BeginStep();
    writer.Put(tm_var, tm_write.data());
    writer.Put(tl_var, tl_write.data());
    writer.Put(cr_var, cr_write.data());
    writer.EndStep();
    writer.Close();
 
    // Old writer to double check results
    // write the events for each processor to their own file
    const fileName exacaPath
    (
        mesh_.time().rootPath()/mesh_.time().globalCaseName()/"ExaCA"
    );
    OFstream os
    (
        exacaPath + "/" + "data_" + Foam::name(Pstream::myProcNo()) + ".csv"
    );

    os << "x,y,z,tm,ts,cr" << endl;

    for(int i=0; i < data.size(); i++)
    {
        int n = data[i].size()-1;

        for(int j=0; j < n; j++)
        {
            os << data[i][j] << ",";
        }
        os << data[i][n] << "\n";
    }
}

bool Foam::functionObjects::ExaCA::end()
{
    //- sort events by cell and in time
    events.shrink();
    sort(events);

    Info<< "Number of solidification events: "
        << returnReduce(events.size(), sumOp<scalar>()) << endl;

    //- map points to cells
    mesh_.time().cpuTimeIncrement();
       
    mapPoints();
    
    Info<< "Successfully mapped points to mesh in: "
        << returnReduce(mesh_.time().cpuTimeIncrement(), maxOp<scalar>()) << " s"
        << endl << endl;

    //- interpolate and write ExaCA data in reduced data format
    const fileName exacaPath
    (
        mesh_.time().rootPath()/mesh_.time().globalCaseName()/"ExaCA"
    );

    mkDir(exacaPath);
    
    interpolate();
    
    Info<< "Successfully interpolated and wrote ExaCA data in: "
        << returnReduce(mesh_.time().cpuTimeIncrement(), maxOp<scalar>()) << " s"
        << endl << endl;

    return true;
}


bool Foam::functionObjects::ExaCA::write()
{
    return true;
}


// ************************************************************************* //
