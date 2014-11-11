BOMD
====
An implementation of the Born-Oppenheimer Molecular Dynamics method based on Hartree-Fock calculations. The quantum mechanical effect of the electrons is included in the calculation of energy and forces for the classical motion of the nuclei.


To compile the following libraries are needed:
- Armadillo
- Libconfig
- HDf5
- Boost

Editing of the code is best done using [QtCreator](http://qt-project.org/downloads).

Example of config file:
```
#-------------------------------
# Example of configuration file
# for the HF program
#-------------------------------
chemicalSystem =
{
    name = "CH4";

    atoms =
    (
        {
        basis = "atom_6_basis_STO-3G.tm";
        position = [0.0, 0.0, 0.0];
        velocity = [0.0, 0.0, 0.0];
        }
        ,
        {
        basis = "atom_1_basis_STO-3G.tm";
        position = [0.        ,  0.        , 2.79679446];
        velocity = [0.0, 0.0, 0.0];
        }
        ,
        {
        basis = "atom_1_basis_STO-3G.tm";
        position = [  0.        ,  2.05980133,  0.         ];
        velocity = [0.0, 0.0, 0.0];
        }
        ,
        {
        basis = "atom_1_basis_STO-3G.tm";
        position = [-1.78384028, -1.02990066,  0.        ];
        velocity = [0.0, 0.0, 0.0];
        }
        ,
        {
        basis = "atom_1_basis_STO-3G.tm";
        position = [ 1.78384028, -1.02990066,  0.       ];
        velocity = [0.0, 0.0, 0.0];
        }
    )

    # Optional (if both are set to zero, there will be equal number of each spin):
    nSpinUpElectrons = 0;
    nSpinDownElectrons = 0;
}


generatorSettings=
{
    dr  = 2.5;

    Nx = 4;
    Ny = 2;
    Nz = 2;

    # "fcc"
    # "c"
    lattice = "c";

    # "none"
    # "uniform"
    # "normal"
    velocityDistribution = "normal";
    temperature = 0.;

};

dynamicSettings =
{
    nSteps = 1000;
    stepSize = 10.;
    boxLength = 20.0;

    #0 = "No BC";
    #1 = "Box";
    #2 = "PeriodicBC"
    BC = 0;
};

modifierSettings =
{
    #0 = "No modifier";
    #1 = "velocityRescaling";

    modifierType = 1;
    velocityRescalingFactor = 0.99;
};


solverSettings =
{
    # 0 = "Restricted HF"
    # 1 = "Unrestricted HF"
    method = 1;

    maxNumOfIteration = 1000000;
    dampingFactor     = 0.5;


    DIISprocedureSettings =
    {
        # 0 = "off"
        # 1 = "on"
        useDIISprocedure = 0;

        iterationLimit = 20;
        nTerms         = 3;

    };

};



fileManagerSettings=
{
    saveResults = 1;
    outputFilePath = "/"
};


analysisSettings=
{
    saveResults = 0;
    outputFilePath = "/"

    saveEnergies = 1;
    dipoleMoment = 1;
    atomicPartialCharge = 1;
    electrostaticPotential = 1;
    chargeDensity = 0;
};

