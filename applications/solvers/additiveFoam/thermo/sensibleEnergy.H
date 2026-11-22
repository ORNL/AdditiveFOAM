#ifndef sensibleEnergy_H
#define sensibleEnergy_H

#include "scalarField.H"
#include "solidLiquidThermo.H"
#include "Polynomial.H"

namespace Foam
{

//- Sensible enthalpy of a cell per volume: the solid-liquid mixture on its
//  path and, below the solidus, the powder the cell may still be. The powder
//  enthalpy is relative to the solid at the solidus, so consolidation there
//  neither absorbs nor releases heat. The enthalpy of the linearized solution
//  is the Taylor step from the last property update.
struct sensibleEnergy
{
    const scalarField& alphaPowder;
    const scalarField& rhoCpCell;
    const scalarField& rhoHCell;
    const solidLiquidThermo& thermo;
    const Polynomial<3>& CpPowder;
    const scalar rho;

    //- Specific heat capacity at T with the solid fraction f of the material
    scalar Cp(const label i, const scalar T, const scalar f) const
    {
        const scalar CpMixture = thermo.Cp(T, f);

        if (T > thermo.solidus())
        {
            return CpMixture;
        }

        return
            (1 - alphaPowder[i])*CpMixture
          + alphaPowder[i]*CpPowder.value(max(T, 300.0));
    }

    //- Specific sensible enthalpy on the path, relative to the solidus
    scalar h(const label i, const scalar T) const
    {
        const scalar Ts = thermo.solidus();
        const scalar hMixture = thermo.h(T);

        if (T > Ts)
        {
            return hMixture;
        }

        const scalar hPowder =
            T >= 300
          ? CpPowder.integral(Ts, T)
          : CpPowder.integral(Ts, 300) + CpPowder.value(300)*(T - 300);

        return (1 - alphaPowder[i])*hMixture + alphaPowder[i]*hPowder;
    }

    scalar rhoHstar(const label i, const scalar T0, const scalar Tstar) const
    {
        return rhoHCell[i] + rhoCpCell[i]*(Tstar - T0);
    }

    scalar rhoH(const label i, const scalar T) const
    {
        return rho*h(i, T);
    }

    scalar rhoCp(const label i, const scalar T) const
    {
        return rho*Cp(i, T, thermo.fs(T));
    }
};

}

#endif
