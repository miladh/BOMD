#include <iostream>
#include <armadillo>
#include <boost/mpi.hpp>

#include <hf.h>
#include <bomd.h>

using namespace arma;
using namespace std;
using namespace hf;
using namespace bomd;


enum solverMethod {
    rhf, uhf
};

enum modifier {
    noModifier, velocityRescaling
};

int main(int argc, char **argv)
{

    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;
    boost::mpi::timer timer;
    timer.restart();

    //read config file---------------------------------------------------------------
    Config cfg;
    for(int p = 0; p < world.size(); p++) {

        world.barrier();
        if(p != world.rank()) {
            continue;
        }
        cfg.readFile("../../../bomd/apps/hydrogen/hydrogenConfig.cfg");
    }
    const Setting & root = cfg.getRoot();


    //Setup system--------------------------------------------------------------------
    string chemicalSystem = root["chemicalSystem"]["name"];
    const Setting &atomsMeta = root["chemicalSystem"]["atoms"];
    double boxLength = root["dynamicSettings"]["boxLength"];
    vector<Atom *> atoms;


    if(atomsMeta.getLength() > 1){
    for(int i = 0; i < atomsMeta.getLength(); i++){
        const Setting &atomMeta = atomsMeta[i];

        string basisFile;
        stringstream basisFilePath;
        rowvec position = zeros<rowvec>(3);
        rowvec velocity = zeros<rowvec>(3);

        atomMeta.lookupValue("basis",basisFile);
        basisFilePath << "infiles/turbomole/"<< basisFile;

        const Setting &pos =  atomMeta["position"];
        const Setting &vel =  atomMeta["velocity"];
        for(int i =0; i < 3; i++){
            position[i] = pos[i];
            velocity[i] = vel[i];
        }

        atoms.push_back(new Atom(basisFilePath.str(), position));
        atoms.at(atoms.size()-1)->setCoreVelocity(velocity);
    }

    }else{
        bomd::Generator generator(&cfg);
        generator.setLattice();
        generator.setVelocity();
        atoms = generator.atoms();
    }
    ElectronicSystem *system = new ElectronicSystem();
    system->addAtoms(atoms);

    //setup solver--------------------------------------------------------------------
    int solverMethod = root["solverSettings"]["method"];
    HFsolver* solver;
    string method;

    switch (solverMethod) {
    case rhf:
        method = "rhf";
        solver = new RHF(system);
        break;

    case uhf:
        method = "uhf";
        solver = new UHF(system);
        break;
    }

    int maxNumOfIteration = root["solverSettings"]["maxNumOfIteration"];
    double dampingFactor = root["solverSettings"]["dampingFactor"];

    solver->setDampingFactor(dampingFactor);
    solver->setMaxNumOfIteration(maxNumOfIteration);

    int useDIISprocdure = root["solverSettings"]["DIISprocedureSettings"]["useDIISprocedure"];
    if(useDIISprocdure){
        int nTerms = root["solverSettings"]["DIISprocedureSettings"]["nTerms"];
        int iterationLimit = root["solverSettings"]["DIISprocedureSettings"]["iterationLimit"];
        solver->useDIISprocedure(nTerms,iterationLimit);
    }


    if(world.rank()==0){
        cout << "---------------------------BOMD------------------------------"  << endl;
        cout << "system:          " << chemicalSystem << endl;
        cout << "method:          " << method << endl;
        cout << "Number of atoms: " << atoms.size() << endl;
    }

    //setup bomd solver and modifiers--------------------------------------------------------------------
    MolecularSystem molecularSystem(&cfg, system, solver);
    int modifierType = root["modifierSettings"]["modifierType"];
    switch (modifierType) {
    case noModifier:
        if(world.rank()==0){
            cout << "No modifiers" <<endl;
        }
        break;

    case velocityRescaling:
        VelocityRescaling* frictionMod;
        frictionMod = new VelocityRescaling(&molecularSystem);
        frictionMod->setRescalingFactor(double(root["modifierSettings"]["velocityRescalingFactor"]));
        molecularSystem.addModifiers(frictionMod);
        if(world.rank()==0){
            cout << " Friction" <<endl;
        }
        break;
    }


    //run solver--------------------------------------------------------------------
    molecularSystem.runDynamics();

    //Save config file-------------------------------------------------------------------------
    if(world.rank() == 0 && int(root["fileManagerSettings"]["saveResults"])){
        string outputFilePath = root["fileManagerSettings"]["outputFilePath"];
        stringstream copyCommand;
        copyCommand << "cp ../../../bomd/apps/hydrogen/hydrogenConfig.cfg" << " " << outputFilePath;
        const char* command = new char[sizeof(copyCommand)];
        command = (copyCommand.str()).c_str();

        int status = std::system( command );
        if(!status){
            cout << "Config file successfully copied!" << endl;
        }else{
            cerr << "Config file cannot be copied!" << endl;
        }
    }


    if(world.rank()==0){
        cout << setprecision(3)
             << "Total elapsed time:  " <<  timer.elapsed()  << "s" << endl;
    }


    return 0;

}


