#ifndef MOLECULARSYSTEM_H
#define MOLECULARSYSTEM_H

#include <iostream>
#include <armadillo>

#include <hf.h>
#include "../fileManager/filemanager.h"

using namespace arma;
using namespace std;
using namespace hf;


namespace bomd{
class MolecularSystem
{
public:
    MolecularSystem(ElectronicSystem *system, HFsolver *solver);
    MolecularSystem(Config *cfg, ElectronicSystem *system, HFsolver *solver);

    void runDynamics();
    void computeForces();

    const mat &energyGradient() const;
    double potentialEnergy() const;

    double frictionConstant() const;
    void setFrictionConstant(double frictionConstant);

    double stepSize() const;
    void setStepSize(double stepSize);

    int nSteps() const;
    void setNSteps(int nSteps);

private:
    Config *m_cfg;
    ElectronicSystem* m_system;
    HFsolver *m_solver;
    GeometricalDerivative* m_GD;
    FileManager* m_outputManager;
    Analyser* m_analyser;
    vector<Atom *> m_atoms;

    int m_nAtoms;
    int m_rank;

    int m_nSteps = 0;
    double m_stepSize = 0.1;
    double m_frictionConstant = 0.0;

    mat m_energyGradient;


    void solveSingleStep();
    void halfKick();
    void systemProperties(int currentTimeStep);
    void freezeAtoms();

    //    void writeLammpsFile(int currentTimeStep);
    void boundaryCheck();
};
}
#endif // MOLECULARSYSTEM_H

