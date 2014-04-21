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
    TEST(energyGradient_H2)
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
        atoms.push_back( new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {-0.5, 0.0, 0.0}));
        atoms.push_back( new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", { 0.5, 0.0, 0.0}));

        ElectronicSystem *system = new ElectronicSystem ();
        system->addAtoms(atoms);

        RHF *solver = new RHF(system);
        BOMD BOSolver(system, solver);

        Atom* atomA = atoms.at(0);
        Atom* atomB = atoms.at(1);


        //Domain
        vec bondLength = linspace(1.0, 4.9, 10);
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

        double h = 1.0E-9;
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
            numericalGradient(x) = -(En - Ep) / (2.0 * h);
        }


        for(uint i = 0; i < numericalGradient.n_elem; i++){
            CHECK_CLOSE(numericalGradient(i), gradient(i), 1e-6);
            //        cout << setprecision(14) << "[" << bondLength(i) << "," << numericalGradient(i) <<"," << gradient(i) << "]," <<endl;

        }

    }
}
