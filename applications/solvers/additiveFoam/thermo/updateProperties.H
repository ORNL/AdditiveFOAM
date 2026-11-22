//- Heat capacity, sensible enthalpy and conductivity of every cell at its
//  state. The powder fraction applies below the solidus only, so
//  consolidation changes no property
#include "thermo/updateEnthalpy.H"
{
    const scalar Ts = solidus.value();

    scalarField& kappaCells = kappa.primitiveFieldRef();

    forAll(T, celli)
    {
        const scalar t = T[celli];
        const scalar f = alphaSolid[celli];

        kappaCells[celli] =
            t > Ts
          ? thermo.kappa(t, f)
          : (1 - alphaPowder[celli])*thermo.kappa(t, f)
          + alphaPowder[celli]*kappaPowder.value(max(t, 300.0));
    }

    kappa.correctBoundaryConditions();
}
