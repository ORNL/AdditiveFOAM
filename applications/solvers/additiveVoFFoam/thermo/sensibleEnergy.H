#ifndef sensibleEnergy_H
#define sensibleEnergy_H

#include "scalarField.H"
#include "solidLiquidThermo.H"
#include "Polynomial.H"

namespace Foam
{

//- Sensible enthalpy of a cell per volume: the metal on its path and the gas,
//  both relative to the solidus. The enthalpy of the linearized solution is
//  the Taylor step from the last property update.
struct sensibleEnergy
{
    const scalarField& alpha;
    const scalarField& rhoCpCell;
    const scalarField& rhoHCell;
    const solidLiquidThermo& thermo;
    const Polynomial<3>& CpGas;
    const scalar rho1;
    const scalar rho2;

    scalar rhoHstar(const label i, const scalar T0, const scalar Tstar) const
    {
        return rhoHCell[i] + rhoCpCell[i]*(Tstar - T0);
    }

    scalar rhoH(const label i, const scalar T) const
    {
        return
            alpha[i]*rho1*thermo.h(T)
          + (1 - alpha[i])*rho2*CpGas.integral(thermo.solidus(), T);
    }

    scalar rhoCp(const label i, const scalar T) const
    {
        return alpha[i]*rho1*thermo.Cp(T) + (1 - alpha[i])*rho2*CpGas.value(T);
    }
};

}

#endif
