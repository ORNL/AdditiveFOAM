//- Face unit interface normal flux, used by the interface compression scheme
{
    const surfaceVectorField gradAlphaf
    (
        fvc::interpolate(fvc::grad(alpha1, "nHat"))
    );

    nHatf = (gradAlphaf/(mag(gradAlphaf) + deltaN)) & mesh.Sf();
}
