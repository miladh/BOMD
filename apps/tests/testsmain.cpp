#include <unittest++/UnitTest++.h>
#include <unittest++/Test.h>
#include <unittest++/TestReporterStdout.h>
#include <unittest++/TestRunner.h>
#include <boost/mpi.hpp>
#include <hf.h>

int main(int argc, char **argv)
{

#if USE_MPI
    boost::mpi::environment env(argc, argv);
#endif
    int result = 0;

 return result;
}
