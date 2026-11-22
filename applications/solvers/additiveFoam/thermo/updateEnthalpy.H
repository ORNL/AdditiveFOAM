//- Heat capacity and sensible enthalpy of every cell at its state. The solid
//  fraction is the state of the cell, which on an isothermal step of the path
//  is not a function of the temperature. primitiveFieldRef() stores the
//  old-time fields before they are overwritten
{
    scalarField& rhoCpCells = rhoCp.primitiveFieldRef();
    scalarField& rhoHCells = rhoH.primitiveFieldRef();

    forAll(T, celli)
    {
        rhoCpCells[celli] = rho.value()*sensible.Cp(celli, T[celli], alphaSolid[celli]);
        rhoHCells[celli] = rho.value()*sensible.h(celli, T[celli]);
    }

    rhoCp.correctBoundaryConditions();
    rhoH.correctBoundaryConditions();
}
