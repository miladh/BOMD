#include <iostream>
#include <armadillo>
#include <boost/mpi.hpp>

#include <hf.h>
#include <bomd.h>

using namespace arma;
using namespace std;

enum solverMethod {
    rhf, uhf
};

hf::ElectronicSystem* setupSystem(string name);
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
    vector<hf::Atom *> atoms;

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

        atoms.push_back(new hf::Atom(basisFilePath.str(), position));
    }

    //    ElectronicSystem *system = new ElectronicSystem();
    //    system->addAtoms(atoms);

    hf::ElectronicSystem *system = setupSystem("H2");

    //setup solver--------------------------------------------------------------------
    int solverMethod = root["solverSettings"]["method"];
    hf::HFsolver* solver;
    string method;

    switch (solverMethod) {
    case rhf:
        method = "rhf";
        solver = new hf::RHF(system);
        break;

    case uhf:
        method = "uhf";
        solver = new hf::UHF(system);
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


    BOMD boSolver(system, solver);
    boSolver.runDynamics();
    double laps = timer.elapsed();


    //Analyzer------------------------------------------------------------------------------
//    Analyser analyser(&cfg, system,solver);
//    analyser.runAnalysis();



    //Save config file-------------------------------------------------------------------------
    if(world.rank() == 0 && int(root["analysisSettings"]["saveResults"])){
        string outputFilePath = root["analysisSettings"]["outputFilePath"];
        stringstream copyCommand;
        copyCommand << "cp ../../../hf/apps/default/defaultConfig.cfg" << " " << outputFilePath;
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
             << "Total elapsed time:  " <<  timer.elapsed()  << "s" << endl
             << " - Computation time: " << laps << "s" << endl
             << " - Analysis time:    " << timer.elapsed() - laps << "s" << endl;
    }



    return 0;

}

hf::ElectronicSystem* setupSystem(string name)
{
    vector<hf::Atom *> atoms;

    if(name =="H2"){
        atoms.push_back(new hf::Atom("infiles/turbomole/atom_1_basis_6-31Gds.tm", { -0.5, 0, 0 }));
        atoms.push_back(new hf::Atom("infiles/turbomole/atom_1_basis_6-31Gds.tm", { 0.2, 0, 0 }));

    }else if(name =="HCl"){
        atoms.push_back(new hf::Atom("infiles/turbomole/atom_17_basis_3-21G.tm", { -0.8, 0, 0 }));
        atoms.push_back(new hf::Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {   0.8, 0, 0 }));

    }else if(name =="Li2"){
        atoms.push_back(new hf::Atom("infiles/turbomole/atom_3_basis_3-21G.tm", {-2.5255, 0.0, 0.0}));
        atoms.push_back(new hf::Atom("infiles/turbomole/atom_3_basis_3-21G.tm", { 2.5255, 0.0, 0.0}));

    }else if(name =="O2"){
        atoms.push_back(new hf::Atom("infiles/turbomole/atom_8_basis_3-21G.tm", {-2.14, 0.0, 0.0}));
        atoms.push_back(new hf::Atom("infiles/turbomole/atom_8_basis_3-21G.tm", { 2.14, 0.0, 0.0}));

    }else if(name =="H2O"){
        atoms.push_back(new hf::Atom("infiles/turbomole/atom_8_basis_3-21G.tm", { 0.0, 0.0, 0.0}));
        atoms.push_back(new hf::Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {2.797, 0.0, 0.0}));
        atoms.push_back(new hf::Atom("infiles/turbomole/atom_1_basis_3-21G.tm", { -1.797*cos((180-104.45) *M_PI/180.0),
                                                                              1.797*sin((180-104.45) *M_PI/180.0), 0.0}));
    }else if(name =="CO2"){
        atoms.push_back(new hf::Atom("infiles/turbomole/atom_8_basis_3-21G.tm", {-2.185, 0.0, 0.0}));
        atoms.push_back(new hf::Atom("infiles/turbomole/atom_8_basis_3-21G.tm", { 2.185, 0.0, 0.0}));
        atoms.push_back(new hf::Atom("infiles/turbomole/atom_6_basis_3-21G.tm", { 0.0, 0.0, 0.0}));

    }else if(name =="CH4"){
        double D = 3.243;
        atoms.push_back(new hf::Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {0.0, 0.0, 0.0}));
        atoms.push_back(new hf::Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {D/sqrt(3), D/sqrt(3), D/sqrt(3)}));
        atoms.push_back(new hf::Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {-D/sqrt(3), -D/sqrt(3), D/sqrt(3)}));
        atoms.push_back(new hf::Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {D/sqrt(3), -D/sqrt(3), -D/sqrt(3)}));
        atoms.push_back(new hf::Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {-D/sqrt(3), D/sqrt(3), -D/sqrt(3)}));

    }else{
        cerr << "unknown system!" << endl;
        exit(0);
    }

    hf::ElectronicSystem *system = new hf::ElectronicSystem();
    system->addAtoms(atoms);

    return system;


}

