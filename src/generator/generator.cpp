#include "generator.h"


bomd::Generator::Generator(Config *cfg):
    m_cfg(cfg),
    m_rank(0)
{
    m_rank = 0;
#if USE_MPI
    boost::mpi::communicator world;
    m_rank = world.rank();
#endif

    const Setting & root = m_cfg->getRoot();
    const Setting &atomMeta = root["chemicalSystem"]["atoms"];
    string basisFile;
    atomMeta[0].lookupValue("basis",basisFile);
    m_basisFilePath << "infiles/turbomole/"<< basisFile;


    m_dr = root["generatorSettings"]["dr"];
    m_Nx = root["generatorSettings"]["Nx"];
    m_Ny = root["generatorSettings"]["Ny"];
    m_Nz = root["generatorSettings"]["Nz"];
}


void bomd::Generator::setLattice()
{
    const Setting & root = m_cfg->getRoot();

    string lattice = root["generatorSettings"]["lattice"];
    if(lattice == "c"){
        cubicLatticeGenerator();
    }else if(lattice == "fcc"){
        fccLatticeGenerator();
    }else{
        cerr << "Unknown lattice!" << endl;
        cerr << lattice << endl;
        exit(0);
    }
}


void bomd::Generator::cubicLatticeGenerator()
{

    rowvec dR = {0., 0., 0.};
    rowvec dr = {m_dr,m_dr, m_dr};

    for (int nZ = 0; nZ < m_Nz; nZ++) {
        dR[2] = nZ * dr[2];
        for (int nY = 0; nY < m_Ny; nY++) {
            dR[1] = nY * dr[1];
            for (int nX = 0; nX < m_Nx; nX++) {
                dR[0] = nX * dr[0];
                m_atoms.push_back(new Atom(m_basisFilePath.str(), dR));
            }
        }
    }
}


void bomd::Generator::fccLatticeGenerator()
{
    rowvec dR = {0., 0., 0.};
    rowvec dr = {m_dr, m_dr, m_dr};
    mat origAtom= zeros<mat>(4,3);


    // FCC atoms in the original unit cell
    origAtom << 0.0 << 0.0 << 0.0 << endr
             << 0.0 << 0.5 << 0.5 << endr
             << 0.5 << 0.0 << 0.5 << endr
             << 0.5 << 0.5 << 0.0 << endr;

    for (int nZ = 0; nZ < m_Nz; nZ++) {
        dR[2] = nZ * dr[2];
        for (int nY = 0; nY < m_Ny; nY++) {
            dR[1] = nY * dr[1];
            for (int nX = 0; nX < m_Nx; nX++) {
                dR[0] = nX * dr[0];
                for (int j=0; j < 4; j++) {
                    m_atoms.push_back(new Atom(m_basisFilePath.str(), dR + dr % origAtom.row(j)));
                }
            }
        }
    }
}


void bomd::Generator::setVelocity()
{

    rowvec sumVelocities = zeros(1,3);
    m_idum = m_idum - time(NULL);
    srand(-m_idum);

    const Setting & root = m_cfg->getRoot();
    string velocityDist = root["generatorSettings"]["velocityDistribution"];


    if(velocityDist == "none"){

    }else if(velocityDist == "uniform"){
        double v = 0.01;
        for(Atom* atom : m_atoms){
            atom->setCoreVelocity({-v+2*v*randu(),-v+2*v*randu(), -v+2*v*randu()});
            sumVelocities +=atom->coreVelocity();
        }
    }else if(velocityDist == "normal"){

        double std = 2.0;
        for(Atom* atom : m_atoms){
            for(int j=0; j < 3 ; j++){
                atom->setCoreVelocity({randn()*std,randn()*std,randn()*std});
            }
            sumVelocities +=atom->coreVelocity();
        }

    }else{
        cerr << "Unknown velocity distribution!" << endl;
        cerr << velocityDist << endl;
        exit(0);
    }


    //Removing initial linear momentum
    for(Atom* atom : m_atoms){
        rowvec v = atom->coreVelocity();
        v -=sumVelocities;
        atom->setCoreVelocity(v);
    }

}
vector<Atom *> bomd::Generator::atoms() const
{
    return m_atoms;
}

