// Implements the tune detector
#ifndef TUNEDETECTOR_H
#define TUNEDETECTOR_H

// MPI Function Wrappers
#include "orbit_mpi.hh"
#include "wrap_mpi_comm.hh"
// ORBIT includes
#include "Bunch.hh"
#include "OrbitConst.hh"
// pyORBIT utils
#include "CppPyWrapper.hh"

// Standard includes
#include <iostream>
#include <vector>
#include <limits>

// Useful macros
#define WARNING_PREFIX "Warning from line " << __LINE__ << " in file " << __FILE__ << ": "

/*****************************************************************************/
// Beam distribution class
class TuneDetector: public OrbitUtils::CppPyWrapper
{
    private:
        int nparticles, trigcnt; // Number of particles and no. of triggers
        std::vector<double> xn, yn, z; // Particle coordinates buffer
        std::vector<int> cntx, cnty, cntz; // Number of positive edges counted

    public:
        TuneDetector(int npart); // Constructor
        ~TuneDetector(); // Destructor
        void reset(int npart); // Reset the counters
        // Track the bunch passing through the monitor
        void trackBunch(Bunch* bunch, double T[][4]);
        // Direct access to the counters
        std::vector<int>& getcntx() {return cntx;};
        std::vector<int>& getcnty() {return cnty;};
        std::vector<int>& getcntz() {return cntz;};
        int getcnttrig() {return trigcnt;};
};

#endif // End of header
