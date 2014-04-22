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

ElectronicSystem* setupSystem(string name);
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
        cfg.readFile("../../../bomd/apps/default/defaultConfig.cfg");
    }
    const Setting & root = cfg.getRoot();


    //Setup system--------------------------------------------------------------------
    string chemicalSystem = root["chemicalSystem"]["name"];
    const Setting &atomsMeta = root["chemicalSystem"]["atoms"];
    vector<Atom *> atoms;

    for(int i = 0; i < atomsMeta.getLength(); i++){
        const Setting &atomMeta = atomsMeta[i];

        string basisFile;
        stringstream basisFilePath;
        rowvec position = zeros<rowvec>(3);

        atomMeta.lookupValue("basis",basisFile);
        basisFilePath << "infiles/turbomole/"<< basisFile;

        const Setting &pos =  atomMeta["position"];
        for(int i =0; i < 3; i++){
            position[i] = pos[i];
        }

        atoms.push_back(new Atom(basisFilePath.str(), position));
    }

    ElectronicSystem *system = new ElectronicSystem();
    system->addAtoms(atoms);

//    ElectronicSystem *system = setupSystem("H4O2");

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



    //run solver--------------------------------------------------------------------
    if(world.rank()==0){
        cout << "---------------------------BOMD------------------------------"  << endl;
        cout << "system:    " << chemicalSystem << endl;
        cout << "method:    " << method << endl;
    }


    BOMD boSolver(&cfg, system, solver);
    boSolver.runDynamics();

    //Save config file-------------------------------------------------------------------------
    if(world.rank() == 0 && int(root["fileManagerSettings"]["saveResults"])){
        string outputFilePath = root["fileManagerSettings"]["outputFilePath"];
        stringstream copyCommand;
        copyCommand << "cp ../../../bomd/apps/default/defaultConfig.cfg" << " " << outputFilePath;
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



ElectronicSystem* setupSystem(string name)
{
    vector<Atom *> atoms;

    if(name =="H2"){
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_6-31Gds.tm", { -0.5, 0, 0 }));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_6-31Gds.tm", { 0.2, 0, 0 }));

    }else if(name =="H4O2"){
        double D  = 7.;
        atoms.push_back(new Atom("infiles/turbomole/atom_8_basis_STO-3G.tm", {0.000, 0.000, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_STO-3G.tm", {2.797, 0.000, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_STO-3G.tm", {-0.45, 1.74, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_8_basis_STO-3G.tm", {D, D, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_STO-3G.tm", {D+2.797, D, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_STO-3G.tm", {-0.45 + D, 1.74 + D, 0.0}));

    }else if(name =="HCl"){
        atoms.push_back(new Atom("infiles/turbomole/atom_17_basis_3-21G.tm", { 0.0, 0, 0 }));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {  2.408644745, 0, 0 }));

    }else if(name =="Li2"){
        atoms.push_back(new Atom("infiles/turbomole/atom_3_basis_3-21G.tm", {-2.5255, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_3_basis_3-21G.tm", { 2.5255, 0.0, 0.0}));

    }else if(name =="O2"){
        atoms.push_back(new Atom("infiles/turbomole/atom_8_basis_3-21G.tm", {-2.14, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_8_basis_3-21G.tm", { 2.14, 0.0, 0.0}));

    }else if(name =="H2O"){
        atoms.push_back(new Atom("infiles/turbomole/atom_8_basis_3-21G.tm", { 0.0, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {2.797, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", { -1.797*cos((180-104.45) *M_PI/180.0),
                                                                              1.7401532877550183, 0.0}));
    }else if(name =="CO2"){
        atoms.push_back(new Atom("infiles/turbomole/atom_8_basis_3-21G.tm", {-2.185, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_8_basis_3-21G.tm", { 2.185, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", { 0.0, 0.0, 0.0}));

    }else if(name =="CH4"){
        double D = 3.243;
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {0.0, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {D/sqrt(3), D/sqrt(3), D/sqrt(3)}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {-D/sqrt(3), -D/sqrt(3), D/sqrt(3)}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {D/sqrt(3), -D/sqrt(3), -D/sqrt(3)}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {-D/sqrt(3), D/sqrt(3), -D/sqrt(3)}));

    }else{
        cerr << "unknown system!" << endl;
        exit(0);
    }

    ElectronicSystem *system = new ElectronicSystem();
    system->addAtoms(atoms);

    return system;


}

