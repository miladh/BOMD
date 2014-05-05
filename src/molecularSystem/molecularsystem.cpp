#include "molecularsystem.h"
#include "../modifier/modifier.h"

using namespace bomd;

MolecularSystem::MolecularSystem(ElectronicSystem *system, HFsolver *solver):
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

MolecularSystem::MolecularSystem(Config *cfg, ElectronicSystem *system, HFsolver *solver):
    MolecularSystem(system, solver)
{
    m_cfg = cfg;
    const Setting & root = m_cfg->getRoot();
    setNSteps(int(root["dynamicSettings"]["nSteps"]));
    setStepSize(double(root["dynamicSettings"]["stepSize"]));

    if(m_rank ==0){
        m_outputManager = new FileManager(m_cfg, m_atoms);
    }

}


void MolecularSystem::runDynamics()
{
    computeForces();
    for(int i = 0; i < m_nSteps; i++){
        if(m_rank == 0 ){
            cout << "MD step:   " << i << endl;
            cout << "-------------------------------------------------------------------------------------"  << endl;

        }
        systemProperties(i);
        solveSingleStep();
    }

    if(m_rank ==0){
        m_outputManager->closeOutput();
    }
}


void MolecularSystem::computeForces()
{
    m_solver->runSolver();
    m_energyGradient = -m_GD->energyGradient();
}

void MolecularSystem::solveSingleStep()
{
    halfKick();
    for(hf::Atom* atom : m_atoms){
        if(!atom->frozen()){
            rowvec corePosition = atom->corePosition() + m_stepSize * atom->coreVelocity();
            atom->setCorePosition(corePosition);
        }
    }
    computeForces();
    halfKick();
    applyModifier();
    boundaryCheck();
}

void MolecularSystem::halfKick()
{
    int i = 0;
    for(hf::Atom* atom : m_atoms){

        if(!atom->frozen()){
        rowvec coreVelocity = atom->coreVelocity() + 0.5 * m_stepSize
                * m_energyGradient.row(i)/(atom->coreMass() * PROTONMASS);

        atom->setCoreVelocity(coreVelocity);
        }
        i++;
    }

}


void MolecularSystem::systemProperties(int currentTimeStep)
{
    int i = currentTimeStep;
    double t = i * m_stepSize;
    double Epot = m_solver->energy();
    double Ekin = 0;
    for (const hf::Atom* atom : m_atoms){
        Ekin += 0.5 * atom->coreMass() * PROTONMASS * dot(atom->coreVelocity(),atom->coreVelocity());
    }

    m_analyser->computeAtomicPartialCharge();
    if(m_rank ==0){
        m_outputManager->writeToFile(i, Ekin, Epot, t);
    }


}


void MolecularSystem::boundaryCheck()
{
    double lim = 10.0;
    for(hf::Atom* atom : m_atoms){
        rowvec corePosition = atom->corePosition();
        rowvec coreVelocity = atom->coreVelocity();

        for(int i = 0; i < 3; i++){
            if( corePosition(i) > lim  || corePosition(i) < -lim){
                coreVelocity(i) *= (-1);
            }
        }

        atom->setCoreVelocity(coreVelocity);
    }
}



void MolecularSystem::addModifiers(Modifier* modifier){
    m_modifiers.push_back(modifier);

}

void MolecularSystem::applyModifier(){

    for(Modifier* modifier: m_modifiers){
        modifier->apply();
    }
}

double MolecularSystem::potentialEnergy() const
{
    return m_solver->energy();
}


const mat& MolecularSystem::energyGradient() const
{
    return m_energyGradient;
}


double MolecularSystem::frictionConstant() const
{
    return m_frictionConstant;
}

void MolecularSystem::setFrictionConstant(double frictionConstant)
{
    m_frictionConstant = frictionConstant;
}

double MolecularSystem::stepSize() const
{
    return m_stepSize;
}

void MolecularSystem::setStepSize(double stepSize)
{
    m_stepSize = stepSize;
}
int MolecularSystem::nSteps() const
{
    return m_nSteps;
}

void MolecularSystem::setNSteps(int nSteps)
{
    m_nSteps = nSteps;
}
vector<Atom *> MolecularSystem::atoms() const
{
    return m_atoms;
}

void MolecularSystem::setAtoms(const vector<Atom *> &atoms)
{
    m_atoms = atoms;
}


