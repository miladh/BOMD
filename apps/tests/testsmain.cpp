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
    int H2Tests = 1;

    UnitTest::TestReporterStdout reporter;
    UnitTest::TestRunner runner(reporter);




    if(H2Tests){
        result += runner.RunTestsIf(UnitTest::Test::GetTestList(), "H2Tests", UnitTest::True(), 0);
    }

    return result;
}



SUITE(H2Tests){
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
        MolecularSystem BOSolver(system, solver);

        Atom* atomA = atoms.at(0);
        Atom* atomB = atoms.at(1);
        rowvec pos;

        BOSolver.computeForces();
//        double dE = atomA->coreAcceleration()(0);
        double E  = BOSolver.potentialEnergy();


        atomB->setCorePosition(X+h);
        double dEn =  atomA->force()(0);
        double En  = BOSolver.potentialEnergy();



        atomB->setCorePosition(X-h);
        BOSolver.computeForces();
        double dEp =  atomA->force()(0);
        double Ep = BOSolver.potentialEnergy();


        atomB->setCorePosition(X+2.0*h);
        BOSolver.computeForces();
        double dEnn = atomA->force()(0);

        atomB->setCorePosition(X-2.0*h);
        BOSolver.computeForces();
        double dEpp =  atomA->force()(0);

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


