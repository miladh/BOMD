#ifndef FILEMANAGER_H
#define FILEMANAGER_H

#include <iostream>
#include <armadillo>
#include <boost/mpi.hpp>
#include <H5Cpp.h>
#include <hf.h>


using namespace arma;
using namespace std;
using namespace H5;
using namespace hf;

namespace bomd{

class FileManager
{
public:
    FileManager(Config *cfg, vector<Atom *> atoms);

    void writeToFile(const int state, const double &kin, const double &pot, const double t);
    void writeLammpsFile(int currentTimeStep);
    void closeOutput();

private:
    Config* m_cfg;
    vector<Atom *> m_atoms;
    int m_nAtoms;
    int m_rank;
    int m_nProcs;
    stringstream m_outputFileName;


    struct AtomAttributes {
        int type;
        char basisType[64];
        double x;
        double y;
        double z;
        double vx;
        double vy;
        double vz;
        int coreCharge;
        double corePartialCharge;
        int frozen;
    };

    H5File *m_output;
    AtomAttributes *m_atomAttributes;
    CompType *m_atomCompound;
    vector <DataSet *>m_dataset;


    void initialize();
};
}
#endif // FILEMANAGER_H
