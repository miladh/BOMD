#ifndef GENERATOR_H
#define GENERATOR_H

#include <armadillo>
#include <iostream>
#include <hf.h>


using namespace arma;
using namespace std;
using namespace hf;

namespace bomd{
class Generator
{
public:
    Generator(Config* cfg, vector<Atom *> *atoms);

    void setLattice();
    void setVelocity();

    vector<Atom *> atoms() const;

    double temperature() const;
    void setTemperature(double temperature);

private:
    Config *m_cfg;
    vector<Atom *>* m_atoms;
    stringstream m_basisFilePath;
    int m_nAtoms;
    int m_rank;
    int m_Nx, m_Ny, m_Nz;
    double m_dr;
    double m_boxLength;
    long m_idum;
    double m_temperature;

    void fccLatticeGenerator();
    void cubicLatticeGenerator();
};
}
#endif // GENERATOR_H
