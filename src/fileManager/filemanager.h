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

namespace bomd{

class FileManager
{
public:
    FileManager(Config *cfg);

    void writeToFile(const int state, vector<hf::Atom *> atoms, const double &kin, const double &pot, const double t);
    void closeOutput();

private:
    int m_rank;
    int m_nProcs;
    Config* m_cfg;
    stringstream m_outputFileName;


    struct AtomAttributes {
        int type;
        char basisType[64];
        double x;
        double y;
        double z;
        int coreCharge;
        double corePartialCharge;
    };

    H5File *m_output;
    AtomAttributes *m_atomAttributes;
    CompType *m_atomCompound;
    vector <DataSet *>m_dataset;


    void initialize();
};
}
#endif // FILEMANAGER_H
