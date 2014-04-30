#include <unittest++/UnitTest++.h>
#include <unittest++/Test.h>
#include <unittest++/TestReporterStdout.h>
#include <unittest++/TestRunner.h>
#include <boost/mpi.hpp>


#include <bomd.h>
#include <hf.h>
#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;
using namespace hf;
using namespace bomd;

int main(int argc, char **argv)
{

#if USE_MPI
    boost::mpi::environment env(argc, argv);
#endif
    int result = 0;

    UnitTest::TestReporterStdout reporter;
    UnitTest::TestRunner runner(reporter);

    result += runner.RunTestsIf(UnitTest::Test::GetTestList(), "DEVELOPMENT", UnitTest::True(), 0);

    return result;
}




SUITE(DEVELOPMENT){

TEST(RHFenergyGradient_H2)
{
    /*
         * test case:   geometrical derivative of energy
         * system:      H2
         * basis:       3-21G
         * method:      RHF
         * source:
         *      numerical differentiation of energy
         * */

    int myRank = 0;
#if USE_MPI
    boost::mpi::communicator world;
    myRank = world.rank();
#endif

    if(myRank == 0){
        cout << "system:    " << "H2" << endl;
        cout << "method:    " << "RHF" << endl;
        cout << "basis:     " << "3-21G" << endl;
    }

    //Initializing the system
    vector<Atom *> atoms;
    atoms.push_back( new Atom("infiles/turbomole/atom_1_basis_STO-3G.tm", {-0.5, 0.0, 0.0}));
    atoms.push_back( new Atom("infiles/turbomole/atom_1_basis_STO-3G.tm", { 0.5, 0.0, 0.0}));

    ElectronicSystem *system = new ElectronicSystem ();
    system->addAtoms(atoms);

    RHF *solver = new RHF(system);
    BOMD BOSolver(system, solver);

    Atom* atomA = atoms.at(0);
    Atom* atomB = atoms.at(1);


    //Domain
    vec bondLength = linspace(1.0, 4.5, 10);
    vec gradient   = 0*bondLength;
    vec numericalGradient = 0*bondLength;

    //Calculate gradient analytically
    for(uint x = 0; x < bondLength.n_elem; x++){
        rowvec X = {bondLength(x) , 0 ,0 };

        atomA->setCorePosition(X * -0.5);
        atomB->setCorePosition(X * 0.5);
        BOSolver.computeForces();

        mat gradE = BOSolver.energyGradient();
        gradient(x) = gradE(1,0);
    }

    double h = 1.0E-7;
    //Calculate gradient numerically
    for(uint x = 0; x < bondLength.n_elem; x++){
        rowvec X = {bondLength(x) , 0 ,0 };
        rowvec dx = {h , 0 ,0 };


        atomA->setCorePosition((X-dx) * -0.5);
        atomB->setCorePosition((X-dx) *  0.5);
        BOSolver.computeForces();
        double Ep = BOSolver.potentialEnergy();

        atomA->setCorePosition((X+dx) * -0.5);
        atomB->setCorePosition((X+dx) *  0.5);
        BOSolver.computeForces();
        double En = BOSolver.potentialEnergy();

        atomA->setCorePosition((X-2.0*dx) * -0.5);
        atomB->setCorePosition((X-2.0*dx) *  0.5);
        BOSolver.computeForces();
        double Epp = BOSolver.potentialEnergy();

        atomA->setCorePosition((X+2.0*dx) * -0.5);
        atomB->setCorePosition((X+2.0*dx) *  0.5);
        BOSolver.computeForces();
        double Enn = BOSolver.potentialEnergy();

        numericalGradient(x) = -(-Enn + 8. * En - 8. * Ep + Epp) / (12.0 * h);
    }


    for(uint i = 0; i < numericalGradient.n_elem; i++){
        CHECK_CLOSE(numericalGradient(i), gradient(i), 1e-8);

    }

}


TEST(UHFenergyGradient_H2)
{
    /*
         * test case:   geometrical derivative of energy
         * system:      H2
         * basis:       3-21G
         * method:      RHF
         * source:
         *      numerical differentiation of energy
         * */

    int myRank = 0;
#if USE_MPI
    boost::mpi::communicator world;
    myRank = world.rank();
#endif

    if(myRank == 0){
        cout << "system:    " << "H2" << endl;
        cout << "method:    " << "RHF" << endl;
        cout << "basis:     " << "3-21G" << endl;
    }

    //Initializing the system
    vector<Atom *> atoms;
    atoms.push_back( new Atom("infiles/turbomole/atom_1_basis_STO-3G.tm", {-0.5, 0.0, 0.0}));
    atoms.push_back( new Atom("infiles/turbomole/atom_1_basis_STO-3G.tm", { 0.5, 0.0, 0.0}));

    ElectronicSystem *system = new ElectronicSystem ();
    system->addAtoms(atoms);

    UHF *solver = new UHF(system);
    solver->setDampingFactor(0.0);
    BOMD BOSolver(system, solver);

    Atom* atomA = atoms.at(0);
    Atom* atomB = atoms.at(1);

    //Domain
    vec bondLength = linspace(1.0, 4.5, 10);
    vec gradient   = 0*bondLength;
    vec numericalGradient = 0*bondLength;

    //Calculate gradient analytically
    for(uint x = 0; x < bondLength.n_elem; x++){
        rowvec X = {bondLength(x) , 0 ,0 };

        atomA->setCorePosition(X * -0.5);
        atomB->setCorePosition(X * 0.5);
        BOSolver.computeForces();

        mat gradE = BOSolver.energyGradient();
        gradient(x) = gradE(1,0);
    }

    double h = 1.0E-7;
    //Calculate gradient numerically
    for(uint x = 0; x < bondLength.n_elem; x++){
        rowvec X = {bondLength(x) , 0 ,0 };
        rowvec dx = {h , 0 ,0 };


        atomA->setCorePosition((X-dx) * -0.5);
        atomB->setCorePosition((X-dx) *  0.5);
        BOSolver.computeForces();
        double Ep = BOSolver.potentialEnergy();

        atomA->setCorePosition((X+dx) * -0.5);
        atomB->setCorePosition((X+dx) *  0.5);
        BOSolver.computeForces();
        double En = BOSolver.potentialEnergy();

        atomA->setCorePosition((X-2.0*dx) * -0.5);
        atomB->setCorePosition((X-2.0*dx) *  0.5);
        BOSolver.computeForces();
        double Epp = BOSolver.potentialEnergy();

        atomA->setCorePosition((X+2.0*dx) * -0.5);
        atomB->setCorePosition((X+2.0*dx) *  0.5);
        BOSolver.computeForces();
        double Enn = BOSolver.potentialEnergy();

        numericalGradient(x) = -(-Enn + 8. * En - 8. * Ep + Epp) / (12.0 * h);
    }


    for(uint i = 0; i < numericalGradient.n_elem; i++){
        CHECK_CLOSE(numericalGradient(i), gradient(i), 1e-8);
    }

}

TEST(vibrationFreq){
    /*
     * test case:   vibration frequency
     * system:      H2
     * basis:       STO-3G
     * method:      RHF
     * source:
     *      numerical differentiation of energy
     * */
    double t0 = 2.418884326505e-17;
    double c= 3e10;
    int myRank = 0;
#if USE_MPI
    boost::mpi::communicator world;
    myRank = world.rank();
#endif

    if(myRank == 0){
        cout << "system:    " << "H2" << endl;
        cout << "method:    " << "RHF" << endl;
        cout << "basis:     " << "STO-3G" << endl;
    }


    //Initializing the system
    rowvec X = {1.380,0.0, 0.0};
    rowvec h = {1e-4, 0.0 , 0.0};
    vector<Atom *> atoms;
    atoms.push_back( new Atom("infiles/turbomole/atom_1_basis_6-31G.tm", {0.0, 0.0, 0.0}));
    atoms.push_back( new Atom("infiles/turbomole/atom_1_basis_6-31G.tm", X));


    ElectronicSystem *system = new ElectronicSystem ();
    system->addAtoms(atoms);

    RHF *solver = new RHF(system);
    BOMD BOSolver(system, solver);

    Atom* atomA = atoms.at(0);
    Atom* atomB = atoms.at(1);
    mat gradE;
    rowvec pos;

    BOSolver.computeForces();
    gradE = BOSolver.energyGradient();
    double dE = gradE(0,0);
    double E  = BOSolver.potentialEnergy();


    atomB->setCorePosition(X+h);
    BOSolver.computeForces();
    gradE = BOSolver.energyGradient();
    double dEn = gradE(0,0);
    double En  = BOSolver.potentialEnergy();



    atomB->setCorePosition(X-h);
    BOSolver.computeForces();
    gradE = BOSolver.energyGradient();
    double dEp = gradE(0,0);
    double Ep = BOSolver.potentialEnergy();


    atomB->setCorePosition(X+2.0*h);
    BOSolver.computeForces();
    gradE = BOSolver.energyGradient();
    double dEnn = gradE(0,0);

    atomB->setCorePosition(X-2.0*h);
    BOSolver.computeForces();
    gradE = BOSolver.energyGradient();
    double dEpp = gradE(0,0);

    double k = (-dEnn + 8. * dEn - 8. * dEp + dEpp) / (12.0 * h(0));
    double ke = (En - 2.* E + Ep)/(h(0)*h(0));


    double mu = PROTONMASS * atomA->coreMass()*atomB->coreMass()/(atomA->coreMass() + atomB->coreMass());

    double freqG = sqrt(k/mu)* 1.0/(2.0*acos((-1)));
    double freqE = sqrt(ke/mu)* 1.0/(2.0*acos((-1)));

    cout << "force constant: " << endl
         << "     grad: "  << k << endl
         << "     enrgy "  << ke << endl;



    cout << "Frequency: " << endl
         << "     grad: "  << freqG/(t0*c)<< endl
         << "     enrgy "  << freqE/(t0*c) << endl;


}
}


