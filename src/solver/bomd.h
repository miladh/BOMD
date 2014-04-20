#ifndef BOMD_H
#define BOMD_H

#include <iostream>
#include <armadillo>

#include <hf.h>

using namespace arma;
using namespace std;

class BOMD
{
public:
    BOMD(hf::ElectronicSystem *system, hf::HFsolver *solver);

    void runDynamics();
    void computeForces();

    const mat &energyGradient() const;
    double potentialEnergy() const;

private:
    hf::ElectronicSystem* m_system;
    hf::HFsolver *m_solver;
    hf::GeometricalDerivative* m_GD;
    vector<hf::Atom *> m_atoms;

    int m_nAtoms;
    int m_nSteps;
    int m_rank;

    double m_dt;
    double m_frictionConstant;

    mat m_energyGradient;

    vec m_time;
    vec m_totalEnergy;
    vec m_kineticEnergy;
    vec m_potentialEnergy;


    void solveSingleStep();
    void initialStep();
    void halfKick();
    void updateCores();

    void writeLammpsFile(int currentTimeStep);
    void writeSystemProperties();


    void systemProperties(int currentTimeStep);
};

#endif // BOMD_H


