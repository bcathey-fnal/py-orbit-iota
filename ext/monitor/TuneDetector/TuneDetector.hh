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
#include <cmath>

// Useful macros
#define WARNING_PREFIX "Warning from line " << __LINE__ << " in file " << __FILE__ << ": "

/*****************************************************************************/
// Simple frequency counter class for single values
class FrequencyCounter
{
    private:
        int sample_cnt, edge_cnt; // Number of samples and positive edges
        double buf; // Buffer to hold previous value
    
    public:
        FrequencyCounter(){reset();}; // Counstructor
        ~FrequencyCounter(){}; // Destructor
        void reset() // Reset the counter
        {
            sample_cnt = edge_cnt = 0; // Counters start at 0
            buf = std::numeric_limits<double>::max(); // Buffer is at positive max
        };
        void sample(double x) // Sample the signal
        {
            sample_cnt++; // Update the number of samples detected
            if(buf < 0.0 && x > 0.0) edge_cnt++; // Number of positive edges
            buf = x; // Save this sample to the buffer
        }
        int getSampleCount(){return sample_cnt;} // Get total number of samples
        int getEdgeCount(){return edge_cnt;}; // Get total number of + edges
};


/*****************************************************************************/
// Beam distribution class
class TuneDetector: public OrbitUtils::CppPyWrapper
{
    private:
        // Number of particles, total number of triggers, number of lattice
        // entrance triggers
        int nparticles, trigcnt_all, trigcnt_first;

        // Frequency counters for all particles in the bunch
        // All triggers of the TuneDetector
        std::vector<FrequencyCounter> fcntx_all, fcnty_all, fcntz_all;
        // Lattice entrance triggers only
        std::vector<FrequencyCounter> fcntx_first, fcnty_first, fcntz_first;
        // General function to calculate tunes
        void gettunes(std::vector<double>& mu,
                std::vector<FrequencyCounter>& fcnt_all,
                std::vector<FrequencyCounter>& fcnt_first, int nturns);

    public:
        TuneDetector(int npart); // Constructor
        ~TuneDetector(); // Destructor
        void reset(int npart); // Reset the counters
        // Track the bunch passing through the monitor
        void trackBunch(Bunch* bunch, double T[][4], bool isfirst=false);

        // Get the particle tunes
        void getTuneX(std::vector<double>& mux, int nturns)
        {
            gettunes(mux, fcntx_all, fcntx_first, nturns);
        };
        void getTuneY(std::vector<double>& muy, int nturns)
        {
            gettunes(muy, fcnty_all, fcnty_first, nturns);
        };
        void getTuneZ(std::vector<double>& muz, int nturns)
        {
            gettunes(muz, fcntz_all, fcntz_first, nturns);
        };
        
        // Direct access to numbers
        int getParticleCount()  {return nparticles;};
        int getTrigCountAll()   {return trigcnt_all;};
        int getTrigCountFirst() {return trigcnt_first;};
};

#endif // End of header
