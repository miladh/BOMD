#include "filemanager.h"

using namespace bomd;
using namespace hf;
FileManager::FileManager(Config *cfg, vector<Atom *> atoms):
    m_cfg(cfg),
    m_atoms(atoms),
    m_nAtoms(atoms.size()),
    m_rank(0),
    m_nProcs(0)
{

#if USE_MPI
    boost::mpi::environment env;
    boost::mpi::communicator world;
    m_rank = world.rank();
    m_nProcs = world.size();
#endif

    const Setting & root = m_cfg->getRoot();
    string output = root["fileManagerSettings"]["outputFilePath"];
    m_outputFileName << output << "/bomdOutput_" << m_rank << ".h5";
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
    m_atomCompound->insertMember("frozen", HOFFSET(AtomAttributes, frozen), PredType::NATIVE_INT);
    //---------------------------------------------------------------------------------------------------------
    initialize();

}

void FileManager::initialize()
{
    const Setting & root = m_cfg->getRoot();
    int    nSteps = root["dynamicSettings"]["nSteps"];
    double stepSize = root["dynamicSettings"]["stepSize"];
    double frictionConstant = double(root["modifierSettings"]["frictionConstant"]);


    Group rootGroup = m_output->openGroup("/");
    Attribute nAtoms_a(rootGroup.createAttribute("nAtoms", PredType::NATIVE_INT, H5S_SCALAR));
    Attribute nSteps_a(rootGroup.createAttribute("nSteps", PredType::NATIVE_INT, H5S_SCALAR));
    Attribute stepSize_a(rootGroup.createAttribute("stepSize", PredType::NATIVE_DOUBLE, H5S_SCALAR));
    Attribute frictionConstant_a(rootGroup.createAttribute("frictionConstant", PredType::NATIVE_DOUBLE, H5S_SCALAR));
    nAtoms_a.write(PredType::NATIVE_INT, &m_nAtoms);
    nSteps_a.write(PredType::NATIVE_INT, &nSteps);
    stepSize_a.write(PredType::NATIVE_DOUBLE, &stepSize);
    frictionConstant_a.write(PredType::NATIVE_DOUBLE, &frictionConstant);


    m_dataset.reserve(nSteps);
    m_atomAttributes = new AtomAttributes[m_nAtoms];

    hsize_t dim[1];
    dim[0] = m_nAtoms;
    DataSpace space(1, dim);

    Group* group = new Group( m_output->createGroup( "/states" ));
    for(int i = 0; i < nSteps; i++){
        stringstream state;
        state << "state" << setw(4) << setfill('0')  << i;
        m_dataset.push_back(new DataSet(group->createDataSet(state.str(), *m_atomCompound, space)));
    }

}

void FileManager::writeToFile(const int state,
                              const double& kin,
                              const double& pot,
                              const double t)
{

    for(int i = 0; i < m_nAtoms; i++) {
        hf::Atom* atom = m_atoms.at(i);
        strcpy(m_atomAttributes[i].basisType, atom->basisType().c_str());
        m_atomAttributes[i].type = atom->atomType();
        m_atomAttributes[i].x = atom->corePosition()(0);
        m_atomAttributes[i].y = atom->corePosition()(1);
        m_atomAttributes[i].z = atom->corePosition()(2);
        m_atomAttributes[i].coreCharge = atom->coreCharge();
        m_atomAttributes[i].corePartialCharge = atom->corePartialCharge();
        m_atomAttributes[i].frozen = int(atom->frozen());
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

    writeLammpsFile(state);


}

void FileManager::closeOutput()
{
    m_output->flush(H5F_SCOPE_GLOBAL);
    m_output->close();
}




void FileManager::writeLammpsFile(int currentTimeStep) {

    const Setting & root = m_cfg->getRoot();
    string output = root["fileManagerSettings"]["outputFilePath"];

    stringstream fileName;
    fileName << output << "/lammps/state" << setw(4) << setfill('0')  << currentTimeStep <<".lmp";
    ofstream lammpsFile(fileName.str(), ios::out | ios::binary);

    // The system boundaries
    double bound = 10.0;
    double xMin = -bound;
    double xMax = bound;
    double yMin = -bound;
    double yMax = bound;
    double zMin = -bound;
    double zMax = bound;
    // Shearing is zero unless the system boundaries are sheared (yes that's "sheared",
    // not "shared")
    double xShear = 0.0;
    double yShear = 0.0;
    double zShear = 0.0;
    // nColumns is the number of data types you want to write. In our case we want to
    // write four - the atom type and the x, y and z components of the position.
    // If you want velocities, forces, etc., just add more columns and write more data.
    int nColumns = 2 + 3 + 3;
    // We could divide the data into chunks by the LAMMPS file format, but we don't - i.e. only
    // use one chunk. The chunk length is then the same as the number of atoms times the number
    // of columns.
    int nChunks = 1;
    int chunkLength = m_nAtoms * nColumns;

    // Write all the above to the lammps file
    lammpsFile.write(reinterpret_cast<const char*>(&currentTimeStep), sizeof(int));
    lammpsFile.write(reinterpret_cast<const char*>(&m_nAtoms), sizeof(int));
    lammpsFile.write(reinterpret_cast<const char*>(&xMin), sizeof(double));
    lammpsFile.write(reinterpret_cast<const char*>(&xMax), sizeof(double));
    lammpsFile.write(reinterpret_cast<const char*>(&yMin), sizeof(double));
    lammpsFile.write(reinterpret_cast<const char*>(&yMax), sizeof(double));
    lammpsFile.write(reinterpret_cast<const char*>(&zMin), sizeof(double));
    lammpsFile.write(reinterpret_cast<const char*>(&zMax), sizeof(double));
    lammpsFile.write(reinterpret_cast<const char*>(&xShear), sizeof(double));
    lammpsFile.write(reinterpret_cast<const char*>(&yShear), sizeof(double));
    lammpsFile.write(reinterpret_cast<const char*>(&zShear), sizeof(double));
    lammpsFile.write(reinterpret_cast<const char*>(&nColumns), sizeof(int));
    lammpsFile.write(reinterpret_cast<const char*>(&nChunks), sizeof(int));
    lammpsFile.write(reinterpret_cast<const char*>(&chunkLength), sizeof(int));

    // Write all the data for each atom to file
    for(int i = 0; i < m_nAtoms; i++) {
        double atomType = m_atoms[i]->atomType();
        double frozen   = m_atoms[i]->frozen();
        rowvec position = m_atoms[i]->corePosition();
        rowvec velocity = m_atoms[i]->coreVelocity();

        lammpsFile.write(reinterpret_cast<const char*>(&atomType), sizeof(double));
        lammpsFile.write(reinterpret_cast<const char*>(&frozen), sizeof(double));
        lammpsFile.write(reinterpret_cast<const char*>(&position(0)), sizeof(double));
        lammpsFile.write(reinterpret_cast<const char*>(&position(1)), sizeof(double));
        lammpsFile.write(reinterpret_cast<const char*>(&position(2)), sizeof(double));
        lammpsFile.write(reinterpret_cast<const char*>(&velocity(0)), sizeof(double));
        lammpsFile.write(reinterpret_cast<const char*>(&velocity(1)), sizeof(double));
        lammpsFile.write(reinterpret_cast<const char*>(&velocity(2)), sizeof(double));

    }
    lammpsFile.close();
}
