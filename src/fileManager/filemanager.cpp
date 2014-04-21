#include "filemanager.h"

using namespace bomd;
FileManager::FileManager(const int nAtoms, const int nSteps, const string& outputFilePath):
    m_rank(0),
    m_nProcs(0)
{

#if USE_MPI
    boost::mpi::environment env;
    boost::mpi::communicator world;
    m_rank = world.rank();
    m_nProcs = world.size();
#endif

    m_outputFileName << outputFilePath << "/bomdOutput_" << m_rank << ".h5";
    m_output = new H5File (m_outputFileName.str(), H5F_ACC_TRUNC);

    //---------------------------------------------------------------------------------------------------------
    m_atomCompound = new CompType(sizeof(AtomAttributes));
    m_atomCompound->insertMember( "type", HOFFSET(AtomAttributes, type), PredType::NATIVE_INT);
    m_atomCompound->insertMember( "basis type", HOFFSET(AtomAttributes, basisType), StrType(PredType::C_S1, 64));
    m_atomCompound->insertMember("x", HOFFSET(AtomAttributes, x), PredType::NATIVE_DOUBLE);
    m_atomCompound->insertMember("y", HOFFSET(AtomAttributes, y), PredType::NATIVE_DOUBLE);
    m_atomCompound->insertMember("z", HOFFSET(AtomAttributes, z), PredType::NATIVE_DOUBLE);
    m_atomCompound->insertMember("core charge", HOFFSET(AtomAttributes, coreCharge), PredType::NATIVE_INT);
    m_atomCompound->insertMember("partial charge", HOFFSET(AtomAttributes, corePartialCharge), PredType::NATIVE_DOUBLE);
   //---------------------------------------------------------------------------------------------------------

   m_dataset.reserve(nSteps);
   m_atomAttributes = new AtomAttributes[nAtoms];

   hsize_t dim[1];
   dim[0] = nAtoms;
   DataSpace space(1, dim);

   Group* group = new Group( m_output->createGroup( "/states" ));
   for(int i = 0; i < nSteps; i++){
       stringstream state;
       state << "state" << setw(4) << setfill('0')  << i;
       m_dataset.push_back(new DataSet(group->createDataSet(state.str(), *m_atomCompound, space)));
   }
}



void FileManager::writeToFile(const int state, vector<hf::Atom *> atoms,
                              const double& kin, const double& pot,
                              const double t)
{

    for(int i = 0; i < signed(atoms.size()); i++) {
        hf::Atom* atom = atoms.at(i);
        strcpy(m_atomAttributes[i].basisType, atom->basisType().c_str());
        m_atomAttributes[i].type = atom->atomType();
        m_atomAttributes[i].x = atom->corePosition()(0);
        m_atomAttributes[i].y = atom->corePosition()(1);
        m_atomAttributes[i].z = atom->corePosition()(2);
        m_atomAttributes[i].coreCharge = atom->coreCharge();
        m_atomAttributes[i].corePartialCharge = atom->corePartialCharge();
    }

    m_dataset[state]->write(m_atomAttributes, *m_atomCompound);

    Attribute kinAtt(m_dataset[state]->createAttribute("kinetic energy", PredType::NATIVE_DOUBLE, H5S_SCALAR));
    Attribute potAtt(m_dataset[state]->createAttribute("potential energy", PredType::NATIVE_DOUBLE, H5S_SCALAR));
    Attribute EtotAtt(m_dataset[state]->createAttribute("total energy", PredType::NATIVE_DOUBLE, H5S_SCALAR));
    Attribute timeAtt(m_dataset[state]->createAttribute("time", PredType::NATIVE_DOUBLE, H5S_SCALAR));

    double Etot = kin + pot;
    kinAtt.write(PredType::NATIVE_DOUBLE, &kin);
    potAtt.write(PredType::NATIVE_DOUBLE, &pot);
    EtotAtt.write(PredType::NATIVE_DOUBLE, &Etot);
    timeAtt.write(PredType::NATIVE_DOUBLE, &t);


}


void FileManager::closeOutput()
{
    m_output->flush(H5F_SCOPE_GLOBAL);
    m_output->close();
}
