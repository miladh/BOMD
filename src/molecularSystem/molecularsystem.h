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
class Modifier;

class MolecularSystem
{
public:
    MolecularSystem(ElectronicSystem *system, HFsolver *solver);
    MolecularSystem(Config *cfg, ElectronicSystem *system, HFsolver *solver);


    void runDynamics();
    void computeForces();

    int nSteps() const;
    double stepSize() const;
    double boxLength() const;
    double potentialEnergy() const;

    vector<Atom *> atoms() const;
    void setAtoms(const vector<Atom *> &atoms);

    void addModifiers(Modifier *modifier);


private:
    Config *m_cfg;
    ElectronicSystem* m_system;
    HFsolver *m_solver;
    GeometricalDerivative* m_GD;
    FileManager* m_outputManager;
    Analyser* m_analyser;
    vector<Atom *> m_atoms;
    vector <Modifier*> m_modifiers;

    int m_nAtoms;
    int m_rank;
    int m_boundaryCondition = 0;

    int m_nSteps = 0;
    int m_MDstep = 0;
    double m_stepSize = 0.1;
    double m_boxLength;

    void solveSingleStep();
    void halfKick();
    void systemProperties(int currentTimeStep);
    void freezeAtoms();
    void boundaryCheck();
    void applyModifier();
    void minimumImageConvention();
};
}
#endif // MOLECULARSYSTEM_H

