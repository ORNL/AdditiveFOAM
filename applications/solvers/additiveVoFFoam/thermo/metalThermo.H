/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                Copyright (C) 2026 Oak Ridge National Laboratory
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
    Foam::metalThermo

Description
    Equilibrium thermodynamics of the metal along the phase-change path:
    solid fraction, sensible enthalpy, heat capacity and conductivity of the
    solid-liquid mixture. Solid properties are evaluated between 300 K and the
    solidus, liquid properties between the liquidus and twice the liquidus.

\*---------------------------------------------------------------------------*/

#ifndef metalThermo_H
#define metalThermo_H

#include "thermoPath.H"
#include "Polynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class metalThermo
{
    // Private Data

        //- Path temperatures, ascending
        scalarField T_;

        //- Solid fraction at the path temperatures
        scalarField fs_;

        //- Integral of the solid fraction from the solidus
        scalarField Ifs_;

        Polynomial<3> CpS_;
        Polynomial<3> CpL_;
        Polynomial<3> kappaS_;
        Polynomial<3> kappaL_;

        scalar Ts_;
        scalar Tl_;

        //- Sensible enthalpy at the liquidus
        scalar hl_;


    // Private Member Functions

        //- Path segment containing T, solidus < T < liquidus
        label segment(const scalar T) const
        {
            label i = 0;
            while (i < T_.size() - 2 && T_[i + 1] < T)
            {
                i++;
            }
            return i;
        }

        //- Integral of the solid fraction from the solidus to T
        scalar integralFs(const scalar T) const
        {
            if (T <= Ts_)
            {
                return 0;
            }
            if (T >= Tl_)
            {
                return Ifs_.last();
            }
            const label i = segment(T);
            return Ifs_[i] + 0.5*(fs_[i] + fs(T))*(T - T_[i]);
        }


public:

    // Constructors

        metalThermo
        (
            const thermoPath& path,
            const dictionary& solidDict,
            const dictionary& liquidDict
        )
        :
            CpS_(solidDict.lookup("Cp")),
            CpL_(liquidDict.lookup("Cp")),
            kappaS_(solidDict.lookup("kappa")),
            kappaL_(liquidDict.lookup("kappa")),
            Ts_(path.solidus()),
            Tl_(path.liquidus())
        {
            labelList order;
            sortedOrder(path.x(), order);

            T_ = scalarField(path.x(), order);
            fs_ = scalarField(path.y(), order);

            for (label i = 1; i < T_.size(); i++)
            {
                if (T_[i] <= T_[i - 1] || fs_[i] > fs_[i - 1])
                {
                    FatalErrorInFunction
                        << "thermoPath must have increasing temperatures and "
                        << "non-increasing solid fractions"
                        << exit(FatalError);
                }
            }

            Ifs_.setSize(T_.size(), 0);
            for (label i = 1; i < T_.size(); i++)
            {
                Ifs_[i] =
                    Ifs_[i - 1] + 0.5*(fs_[i - 1] + fs_[i])*(T_[i] - T_[i - 1]);
            }

            hl_ =
                CpL_.value(Tl_)*(Tl_ - Ts_)
              + (CpS_.value(Ts_) - CpL_.value(Tl_))*Ifs_.last();
        }


    // Member Functions

        scalar solidus() const
        {
            return Ts_;
        }

        scalar liquidus() const
        {
            return Tl_;
        }

        //- Equilibrium solid fraction
        scalar fs(const scalar T) const
        {
            if (T <= Ts_)
            {
                return 1;
            }
            if (T >= Tl_)
            {
                return 0;
            }
            const label i = segment(T);
            return
                fs_[i]
              + (fs_[i + 1] - fs_[i])*(T - T_[i])/(T_[i + 1] - T_[i]);
        }

        //- Derivative of the solid fraction with temperature
        scalar dfsdT(const scalar T) const
        {
            if (T <= Ts_ || T >= Tl_)
            {
                return 0;
            }
            const label i = segment(T);
            return (fs_[i + 1] - fs_[i])/(T_[i + 1] - T_[i]);
        }

        //- Specific heat capacity
        scalar Cp(const scalar T) const
        {
            if (T <= Ts_)
            {
                return CpS_.value(max(T, 300.0));
            }
            if (T >= Tl_)
            {
                return CpL_.value(min(T, 2*Tl_));
            }
            const scalar f = fs(T);
            return f*CpS_.value(Ts_) + (1 - f)*CpL_.value(Tl_);
        }

        //- Sensible enthalpy relative to the solidus
        scalar h(const scalar T) const
        {
            if (T <= Ts_)
            {
                return
                    T >= 300
                  ? CpS_.integral(Ts_, T)
                  : CpS_.integral(Ts_, 300) + CpS_.value(300)*(T - 300);
            }
            if (T >= Tl_)
            {
                return
                    hl_
                  + (
                        T <= 2*Tl_
                      ? CpL_.integral(Tl_, T)
                      : CpL_.integral(Tl_, 2*Tl_)
                      + CpL_.value(2*Tl_)*(T - 2*Tl_)
                    );
            }
            return
                CpL_.value(Tl_)*(T - Ts_)
              + (CpS_.value(Ts_) - CpL_.value(Tl_))*integralFs(T);
        }

        //- Thermal conductivity
        scalar kappa(const scalar T) const
        {
            const scalar f = fs(T);
            return
                f*kappaS_.value(min(max(T, 300.0), Ts_))
              + (1 - f)*kappaL_.value(min(max(T, Tl_), 2*Tl_));
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
