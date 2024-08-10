// Implements the tune detector
#ifndef TUNEDETECTOR_H
#define TUNEDETECTOR_H

// ORBIT includes
#include "Bunch.hh"
#include "OrbitConst.hh"
// pyORBIT utils
#include "CppPyWrapper.hh"

// Standard includes
#include <iostream>
#include <vector>
#include <limits>
#include <numeric>
#include <algorithm>
#include <cmath>

// FFTW library header
#include "fftw3.h"

// Useful macros
#define WARNING_PREFIX "Warning from line " << __LINE__ << " in file " << __FILE__ << ": "

/*****************************************************************************/
// 1D FFT storage
class FFT1D
{
    private:
        int bufsize;
        fftw_complex *fft_in, *fft_out; // Declare the data buffers for DFT
        fftw_plan fft_p; // Declare a plan to compute the DFT
        void clearfftw(); // Clear fftw variables

    public:
        FFT1D(int nsamples); // Constructor
        ~FFT1D(); // Destructor
        void reset(int nsamples); // Reset the calculator
        double getNormalizedFrequency(std::vector<double> &buffer);
};

/*****************************************************************************/
// Frequency detector class for single values
class FrequencyDetector
{
    private:
        // Size of buffer, number of samples and positive edges
        int bufsize, sample_cnt, edge_cnt;
        double prev; // Hold previous value
        std::vector<double> buffer; // Buffer to hold all the values
        FFT1D *calc_fftw_dft_1d; // Pointer to a 1D FFT calculator
    
    public:
        FrequencyDetector(int nsamples=0, FFT1D *dft_1d=NULL); // Counstructor
        ~FrequencyDetector(){}; // Destructor
        void reset(int nsamples=0, FFT1D *dft_1d=NULL); // Reset
        void sample(double x); // Sample the signal
        int getSampleCount(){return sample_cnt;} // Get total number of samples
        int getEdgeCount(){return edge_cnt;}; // Get total number of + edges
        int getBufferSize(){return bufsize;} // Get the buffer size
        double getNormalizedFrequency(); // Get the normalized frequency
};


/*****************************************************************************/
// Beam distribution class
class TuneDetector: public OrbitUtils::CppPyWrapper
{
    private:
        // Number of particles, total number of triggers,
        // number of lattice and entrance triggers
        int nparticles, trigcnt_all, trigcnt_first;

        FFT1D *calc_fftw_dft_1d; // Pointer to a 1D FFT calculator

        // Frequency counters for all particles in the bunch
        // All triggers of the TuneDetector
        std::vector<FrequencyDetector> fcntx_all, fcnty_all, fcntz_all;
        // Lattice entrance triggers only
        std::vector<FrequencyDetector> fcntx_first, fcnty_first, fcntz_first;
        // General function to calculate tunes
        void gettunes(std::vector<double>& mu,
                std::vector<FrequencyDetector>& fcnt_all,
                std::vector<FrequencyDetector>& fcnt_first,
                bool usefirst);

    public:
        TuneDetector(int npart, int nturns); // Constructor
        ~TuneDetector(){delete calc_fftw_dft_1d;}; // Destructor
        void reset(int npart, int nturns); // Reset the counters
        // Track the bunch passing through the monitor
        void trackBunch(Bunch* bunch, double T[][5], bool isfirst=false);

        // Get the particle tunes
        void getTuneX(std::vector<double>& mux, bool usefirst=false)
        {
            gettunes(mux, fcntx_all, fcntx_first, usefirst);
        };
        void getTuneY(std::vector<double>& muy, bool usefirst=false)
        {
            gettunes(muy, fcnty_all, fcnty_first, usefirst);
        };
        void getTuneZ(std::vector<double>& muz, bool usefirst=false)
        {
            gettunes(muz, fcntz_all, fcntz_first, usefirst);
        };
        
        // Direct access to numbers
        int getParticleCount()  {return nparticles;};
        int getTrigCountAll()   {return trigcnt_all;};
        int getTrigCountFirst() {return trigcnt_first;};
};

#endif // End of header
