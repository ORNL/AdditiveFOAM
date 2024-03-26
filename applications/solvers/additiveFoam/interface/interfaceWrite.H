const fileName interfacePath
(
    runTime.rootPath()/runTime.globalCaseName()/"solidificationData"
);

mkDir(interfacePath);

OFstream os
(
    interfacePath + "/" + "data_" + Foam::name(Pstream::myProcNo()) + ".csv"
);

for(int i=0; i < interfaceData.size(); i++)
{
    int n = interfaceData[i].size()-1;

    for(int j=0; j < n; j++)
    {
        os << interfaceData[i][j] << ",";
    }
    os << interfaceData[i][n] << endl;
}
