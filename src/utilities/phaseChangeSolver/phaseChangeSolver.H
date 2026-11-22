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

Class
    Foam::phaseChangeSolver

Description
    Phase-change iteration of the energy equation on a thermoPath.

    Within a time step the latent heat is linearized about the current
    temperature (Newton), the solver's energy equation is solved with the
    linearized latent source, and every cell is then put back on the path at
    the enthalpy of the linearized solution. The iteration ends when the
    largest change of the solid fraction in a cell falls below the tolerance.

    The linearized latent heat is written with the time-derivative operator
    of the energy equation, fvm::ddt(dFdT, T) - fvc::ddt(dFdT, T), so the
    iteration is consistent with whichever ddt scheme the temperature uses;
    a solver that integrates explicitly requests the Euler operator instead.

    setFromEnthalpy() puts every cell on the path at a given enthalpy per
    volume, which keeps the enthalpy of a cell authoritative across a mesh
    change: the solver stores it before the change, the mesh maps it, and the
    temperature and solid fraction are recovered from it afterwards.

    The solid volume fraction of a cell is alpha*F(T), where F is the
    thermoPath and alpha an optional volume fraction of the material that
    changes phase (alpha.metal in additiveVoFFoam), one everywhere when no
    field is given. The solver keeps its own form of the sensible energy
    equation and supplies the sensible enthalpy of a cell to project() as an
    object with three member functions, all per unit volume:

        rhoHstar(celli, T0, Tstar)  sensible enthalpy of the linearized
                                    solution Tstar about T0
        rhoH(celli, T)             sensible enthalpy at T
        rhoCp(celli, T)            its derivative with temperature

Usage
    \verbatim
        phaseChangeSolver phaseChange(T, alphaSolid, path, rho*Lf, &alpha);

        phaseChange.begin(nCorrectors, tolerance);

        while (phaseChange.loop())
        {
            fvScalarMatrix TEqn
            (
                sensibleTerms == phaseChange.latentSource()
            );

            TEqn.solve();

            phaseChange.project(sensible);
        }
    \endverbatim

SourceFiles
    phaseChangeSolver.C
    phaseChangeSolverTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef phaseChangeSolver_H
#define phaseChangeSolver_H

#include "volFields.H"
#include "fvMatrices.H"
#include "thermoPath.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class phaseChangeSolver Declaration
\*---------------------------------------------------------------------------*/

class phaseChangeSolver
{
    // Private Data

        const fvMesh& mesh_;

        //- Temperature
        volScalarField& T_;

        //- Solid volume fraction of the cell
        volScalarField& alphaSolid_;

        //- Phase-change path: solid fraction of the material against
        //  temperature
        const thermoPath& path_;

        //- Latent heat per volume of the material
        const dimensionedScalar rhoLf_;

        //- Optional volume fraction of the material; one everywhere when null
        const volScalarField* alphaPtr_;

        //- Temperature of the linearization
        volScalarField T0_;

        //- Derivative of the solid fraction with temperature at T0
        volScalarField dFdT_;

        //- Solid fraction before the solve
        volScalarField alphaSolid0_;

        //- Maximum number of iterations in the step
        label nCorr_;

        //- Tolerance on the change of the solid fraction in a cell
        scalar tol_;

        //- Iterations taken in the step
        label iter_;

        //- Largest change of the solid fraction in a cell by the last solve
        scalar residual_;

        bool converged_;


    // Private Member Functions

        //- Linearize the solid fraction about the current temperature
        void linearize();

        //- Volume fraction of the material in a cell
        inline scalar alpha(const label celli) const
        {
            return alphaPtr_ ? (*alphaPtr_)[celli] : 1;
        }

        //- Put one cell on the path at the enthalpy rhoHstar, starting
        //  from TCell with the root bracketed from [TLo, THi]; returns the
        //  iterations taken, maxCellIter when the cap was reached
        template<class SensibleEnergy>
        label solveCell
        (
            const label celli,
            const scalar rhoHstar,
            const scalar TLo,
            const scalar THi,
            const SensibleEnergy& sensible,
            scalar& TCell,
            scalar& alphaSolid
        ) const;


public:

    // Static Data

        //- Cap on the iterations that put one cell back on the path
        static const label maxCellIter = 60;


    // Constructors

        phaseChangeSolver
        (
            volScalarField& T,
            volScalarField& alphaSolid,
            const thermoPath& path,
            const dimensionedScalar& rhoLf,
            const volScalarField* alpha = nullptr
        );

        //- Disallow default bitwise copy construction
        phaseChangeSolver(const phaseChangeSolver&) = delete;


    // Member Functions

        //- Start the iteration of a time step
        void begin(const label nCorr, const scalar tol);

        //- Linearize for the next solve; false when the step is converged or
        //  the iterations are exhausted, with a warning in the latter case
        bool loop();

        //- Latent heat terms of the energy equation for the linearization,
        //  rhoLf*(ddt(alphaSolid) + fvm::ddt(dFdT, T) - fvc::ddt(dFdT, T)),
        //  with the ddt scheme of T, or the Euler scheme when requested
        tmp<fvScalarMatrix> latentSource(const bool eulerDdt = false) const;

        //- Put every cell on the path at the enthalpy of the solved
        //  linearization, and update the residual
        template<class SensibleEnergy>
        void project(const SensibleEnergy& sensible);

        //- Set the temperature and solid fraction of every cell from its
        //  enthalpy per volume, for instance after a mesh change has mapped
        //  the enthalpy
        template<class SensibleEnergy>
        void setFromEnthalpy
        (
            const volScalarField& rhoH,
            const SensibleEnergy& sensible
        );

        //- Iterations taken in the step
        inline label nIter() const
        {
            return iter_;
        }

        inline scalar residual() const
        {
            return residual_;
        }

        inline bool converged() const
        {
            return converged_;
        }

        //- Temperature of the linearization
        inline const volScalarField& T0() const
        {
            return T0_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const phaseChangeSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "phaseChangeSolverTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
