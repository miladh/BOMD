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
    m_boxLength =

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

    m_stepSize   = root["dynamicSettings"]["stepSize"];
    m_nSteps     = root["dynamicSettings"]["nSteps"];
    m_boxLength  = root["dynamicSettings"]["boxLength"];
    m_boundaryCondition = root["dynamicSettings"]["BC"];


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
        m_MDstep = i;
    }

    if(m_rank ==0){
        m_outputManager->closeOutput();
    }
}


void MolecularSystem::computeForces()
{
    if(m_boundaryCondition == 2){
        minimumImageConvention();
    }
    else{
        m_solver->runSolver();
        mat energyGradient = -m_GD->energyGradient();
        for(int i = 0; i < m_nAtoms; i++){
            m_atoms.at(i)->addForce(energyGradient.row(i));
        }
    }
}


void MolecularSystem::minimumImageConvention()
{
    mat m_corePostions = zeros(m_nAtoms, 3);
    for(int i = 0; i < m_nAtoms; i++){
        m_corePostions.row(i) = m_atoms.at(i)->corePosition();
    }


    cube dR = zeros(m_nAtoms, 3,  m_nAtoms);
    for(int i = 0; i < m_nAtoms; i++){
        for(int j = 0; j < m_nAtoms; j++){
            const rowvec& dr = m_atoms.at(j)->corePosition() - m_atoms.at(i)->corePosition();

            for(int k = 0; k < 3; k++){
                dR.slice(i)(j,k) = dr(k)  - m_boxLength * round(dr(k) /m_boxLength);
            }
        }
    }

    for(int i = 0; i < m_nAtoms; i++){
        for(int j = 0; j < m_nAtoms; j++){
            rowvec pos =  m_corePostions.row(i) + dR.slice(i).row(j);
            m_atoms.at(j)->setCorePosition(pos);
        }

        m_solver->runSolver();
        m_atoms.at(i)->addForce(-m_GD->energyGradient().row(i));
    }

    for(int k = 0; k < m_nAtoms; k++){
        m_atoms.at(k)->setCorePosition(m_corePostions.row(k));
    }

}


void MolecularSystem::boundaryCheck()
{

    if(m_boundaryCondition == 0);

    else if (m_boundaryCondition == 1){
        for(hf::Atom* atom : m_atoms){
            rowvec corePosition = atom->corePosition();
            rowvec coreVelocity = atom->coreVelocity();

            for(int i = 0; i < 3; i++){
                if( corePosition(i) > m_boxLength  || corePosition(i) < 0){
                    coreVelocity(i) *= (-1);
                }
            }
            atom->setCoreVelocity(coreVelocity);
        }
    }
    else if(m_boundaryCondition == 2){
        for(hf::Atom* atom : m_atoms){
            rowvec corePosition = atom->corePosition();
            for(int i = 0; i < 3; i++){
                corePosition(i) = fmod((corePosition(i) + 2 * m_boxLength), m_boxLength);
            }

            atom->setCorePosition(corePosition);
        }
    }

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
    boundaryCheck();
    computeForces();
    halfKick();
    applyModifier();
}

void MolecularSystem::halfKick()
{
    int i = 0;
    for(hf::Atom* atom : m_atoms){

        if(!atom->frozen()){
            rowvec coreVelocity = atom->coreVelocity() + 0.5 * m_stepSize
                    * atom->force()/(atom->coreMass() * PROTONMASS);

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




void MolecularSystem::addModifiers(Modifier* modifier){
    m_modifiers.push_back(modifier);

}

double MolecularSystem::boxLength() const
{
    return m_boxLength;
}

void MolecularSystem::applyModifier()
{
    for(Modifier* modifier: m_modifiers){
        modifier->apply();
    }
}

double MolecularSystem::potentialEnergy() const
{
    return m_solver->energy();
}

double MolecularSystem::stepSize() const
{
    return m_stepSize;
}

int MolecularSystem::nSteps() const
{
    return m_nSteps;
}

vector<Atom *> MolecularSystem::atoms() const
{
    return m_atoms;
}

void MolecularSystem::setAtoms(const vector<Atom *> &atoms)
{
    m_atoms = atoms;
}


