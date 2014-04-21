#include "bomd.h"

using namespace bomd;

BOMD::BOMD(ElectronicSystem *system, HFsolver *solver):
    m_system(system),
    m_solver(solver),
    m_atoms(system->atoms()),
    m_nAtoms(system->nAtoms())
{
    m_GD = new hf::GeometricalDerivative(m_system, m_solver);
    m_analyser = new Analyser(m_system, m_solver);

    m_rank = 0;
#if USE_MPI
    boost::mpi::communicator world;
    m_rank = world.rank();
#endif

}

BOMD::BOMD(Config *cfg, ElectronicSystem *system, HFsolver *solver):
    BOMD(system, solver)
{
    m_cfg = cfg;
    const Setting & root = m_cfg->getRoot();
    setNSteps(int(root["dynamicSettings"]["nSteps"]));
    setStepSize(double(root["dynamicSettings"]["stepSize"]));
    setFrictionConstant(double(root["dynamicSettings"]["firctionConstant"]));

    string outputFilePath = root["analysisSettings"]["outputFilePath"];
    if(m_rank ==0){
        m_outputManager = new FileManager(m_system->nAtoms(),m_nSteps, outputFilePath);
    }
}



void BOMD::runDynamics()
{
    computeForces();
    for(int i = 0; i < m_nSteps; i++){
        if(m_rank == 0 ){
            cout << "MD step:   " << i << endl;
            cout << "-------------------------------------------------------------------------------------"  << endl;

        }
        solveSingleStep();
        writeLammpsFile(i);
        systemProperties(i);
    }

    if(m_rank ==0){
        m_outputManager->closeOutput();
    }
}


void BOMD::computeForces()
{
    m_solver->runSolver();
    m_energyGradient = -m_GD->energyGradient();
}

void BOMD::solveSingleStep()
{
    halfKick();
    for(hf::Atom* atom : m_atoms){
        rowvec corePosition = atom->corePosition() + m_stepSize * atom->coreVelocity();
        atom->setCorePosition(corePosition);
    }
    computeForces();
    halfKick();
}

void BOMD::halfKick()
{
    int i = 0;
    for(hf::Atom* atom : m_atoms){
        rowvec coreVelocity = atom->coreVelocity() + 0.5 * m_stepSize
                            * (m_energyGradient.row(i)/atom->coreMass()
                            -  m_frictionConstant * atom->coreVelocity());

        atom->setCoreVelocity(coreVelocity);
        i++;
    }

}


void BOMD::systemProperties(int currentTimeStep)
{
    int i = currentTimeStep;
    double t = i * m_stepSize;
    double Epot = m_solver->energy();;
    double Ekin = 0;
    for (const hf::Atom* atom : m_atoms){
        Ekin += 0.5 * atom->coreMass() * dot(atom->coreVelocity(),atom->coreVelocity());
    }

    m_analyser->computeAtomicPartialCharge();
    if(m_rank ==0){
    m_outputManager->writeToFile(i,m_atoms, Ekin, Epot, t);
    }


}


double BOMD::potentialEnergy() const
{
    return m_solver->energy();
}


const mat& BOMD::energyGradient() const
{
    return m_energyGradient;
}


double BOMD::frictionConstant() const
{
    return m_frictionConstant;
}

void BOMD::setFrictionConstant(double frictionConstant)
{
    m_frictionConstant = frictionConstant;
}

double BOMD::stepSize() const
{
    return m_stepSize;
}

void BOMD::setStepSize(double stepSize)
{
    m_stepSize = stepSize;
}
int BOMD::nSteps() const
{
    return m_nSteps;
}

void BOMD::setNSteps(int nSteps)
{
    m_nSteps = nSteps;
}


void BOMD::writeLammpsFile(int currentTimeStep) {

    stringstream outStepName;
    outStepName <<"/home/milad/kurs/qmd/state" << setw(4) << setfill('0')  << currentTimeStep <<".lmp";
    ofstream lammpsFile(outStepName.str(), ios::out | ios::binary);

    // The system boundaries
    double xMin = -5.0;
    double xMax = 5.0;
    double yMin = -5.0;
    double yMax = 5.0;
    double zMin = -5.0;
    double zMax = 5.0;
    // Shearing is zero unless the system boundaries are sheared (yes that's "sheared",
    // not "shared")
    double xShear = 0.0;
    double yShear = 0.0;
    double zShear = 0.0;
    // nColumns is the number of data types you want to write. In our case we want to
    // write four - the atom type and the x, y and z components of the position.
    // If you want velocities, forces, etc., just add more columns and write more data.
    int nColumns = 1 + 3 + 3 + 3;
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
        rowvec position = m_atoms[i]->corePosition();
        rowvec velocity = m_atoms[i]->coreVelocity();

        lammpsFile.write(reinterpret_cast<const char*>(&atomType), sizeof(double));
        lammpsFile.write(reinterpret_cast<const char*>(&position(0)), sizeof(double));
        lammpsFile.write(reinterpret_cast<const char*>(&position(1)), sizeof(double));
        lammpsFile.write(reinterpret_cast<const char*>(&position(2)), sizeof(double));
        lammpsFile.write(reinterpret_cast<const char*>(&velocity(0)), sizeof(double));
        lammpsFile.write(reinterpret_cast<const char*>(&velocity(1)), sizeof(double));
        lammpsFile.write(reinterpret_cast<const char*>(&velocity(2)), sizeof(double));
        lammpsFile.write(reinterpret_cast<const char*>(&m_energyGradient(i,0)), sizeof(double));
        lammpsFile.write(reinterpret_cast<const char*>(&m_energyGradient(i,1)), sizeof(double));
        lammpsFile.write(reinterpret_cast<const char*>(&m_energyGradient(i,2)), sizeof(double));

    }
    lammpsFile.close();
}

