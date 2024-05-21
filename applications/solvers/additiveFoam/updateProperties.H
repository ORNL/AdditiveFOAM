forAll(mesh.cells(), cellI)
{
    // update phase fractions
    scalar alpha1_ = alpha1[cellI];
    scalar alpha2_ = 1.0 - alpha1[cellI];
    
    if (alpha1_ < 1.0)
    {
    	alpha3[cellI] = 0.0;
    }
    scalar alpha3_ = alpha3[cellI];

    // update temperature ranges for phase properties
    scalar T1 = min(max(T[cellI], 300.0), Tsol.value());
    scalar T2 = min(max(T[cellI], Tliq.value()), 2.0*Tliq.value());

    // update solid-liquid-powder mixture properties    
    kappa[cellI] = 
        (1.0 - alpha3_)*(alpha1_*kappa1.value(T1) + alpha2_*kappa2.value(T2))
      + alpha3_*kappa3.value(T1);
      
    Cp[cellI] =
    	(1.0 - alpha3_)*(alpha1_*Cp1.value(T1) + alpha2_*Cp2.value(T2))
      + alpha3_*Cp3.value(T1);
}

alpha3.correctBoundaryConditions();
Cp.correctBoundaryConditions();
kappa.correctBoundaryConditions();
