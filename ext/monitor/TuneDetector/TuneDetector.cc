// Implements the tune detector

#include "TuneDetector.hh"

using namespace OrbitUtils;

/*****************************************************************************/
// Implementation of the FFT1D class
// Constructor
FFT1D::FFT1D(int nsamples)
{
    fft_in = fft_out = NULL;
    fft_p = NULL;
    bufsize = 0;
    reset(nsamples);
}

// Destructor
FFT1D::~FFT1D()
{
    clearfftw(); // Clear FFTW variables
}

// Clear FFTW variables
void FFT1D::clearfftw()
{
    // Free the FFTW variables
    if(fft_p)fftw_destroy_plan(fft_p);
    if(fft_in)fftw_free(fft_in);
    if(fft_out)fftw_free(fft_out);
}

void FFT1D::reset(int nsamples)
{
    if(nsamples == bufsize) return; // Should we reset the calculator?
    bufsize = nsamples; // Save the number of samples
    clearfftw(); // Clear FFTW variables
    // Allocate space for holding DFT data
    fft_in = fftw_alloc_complex(nsamples);
    fft_out = fftw_alloc_complex(nsamples);
    // Plan how to compute the DFT
    fft_p = fftw_plan_dft_1d(nsamples, fft_in, fft_out,
                             FFTW_FORWARD, FFTW_MEASURE);
}

// Get the normalized frequency for a gven signal
double FFT1D::getNormalizedFrequency(std::vector<double> &buffer)
{
    int i, maxindex; // Iteration index, index of DFT with maximum amplitude
    // DC offset, amplitude of DFT and its maximum
    double dc_offset, ampsqr, maxampsqr = 0.0;
    if(bufsize != buffer.size())
        std::cout << WARNING_PREFIX << "Signal buffer size must equal allocated FFTW buffer size." << std::endl;

    // Find the DC offset for the signal in the buffer which is the mean of
    // of the signal.
    dc_offset = std::accumulate(buffer.begin(), buffer.end(), 0.0)/bufsize;
    // Copy the data into the FFTW buffer
    for(i=0; i<bufsize; i++)
    {
        fft_in[i][0] = buffer[i] - dc_offset;
        fft_in[i][1] = 0.0;
    }
    fftw_execute(fft_p); // Compute the DFT using FFTW
    // Analyze the DFT to find the dominant frequency
    // Get rid of synchrotron ocillation somehow??
    for(i=50; i<bufsize/2; i++)
    {
        // Calculate the DFT amplitude
        ampsqr = fft_out[i][0]*fft_out[i][0] + fft_out[i][1]*fft_out[i][1];
        if(maxampsqr < ampsqr) // Check if its larger than the saved value
        {
            maxindex = i;
            maxampsqr = ampsqr;
        }
    }

    /*if(maxindex < 100)
    {
        std::cout << WARNING_PREFIX << "Data dump follows:" << std::endl;
        for(i=0; i<bufsize; i++)
            std::cout << buffer[i] << "\t" << fft_in[i][0] << "\t" << 
                fft_out[i][0]*fft_out[i][0] + fft_out[i][1]*fft_out[i][1] << std::endl;
        std::cout << std::endl;
    }*/

    return (double)maxindex/bufsize; // Return the normalized frequency
}


/*****************************************************************************/
// Implementation of the FrequencyDetector class
// Constructor
FrequencyDetector::FrequencyDetector(int nsamples, FFT1D *dft_1d)
{
    // Save the instance of the 1D DFT convenience class
    calc_fftw_dft_1d = dft_1d;
    bufsize = 0; // By default, the buffer size is 0
    reset(nsamples); // Simply call reset!
}

// Reset the counters and the buffers
void FrequencyDetector::reset(int nsamples, FFT1D *dft_1d)
{
    int i; // Iteration index
    // Reset the simple frequency couner
    sample_cnt = edge_cnt = 0; // Counters start at 0
    // Assume initial value is at positive max
    prev = std::numeric_limits<double>::max();

    // Save the instance of the 1D DFT convenience class
    calc_fftw_dft_1d = dft_1d;
    
    if(nsamples == 0)
    {
        bufsize = 0;
    }
    else if(nsamples != bufsize)
    {
        buffer.resize(nsamples); // Resize the data buffer
        bufsize = nsamples; // Number of samples in the buffer
        // Initialize the buffer to 0.0
        for(i=0; i<bufsize; i++)buffer[i] = 0.0;
    }
} 

// Sample the signal
void FrequencyDetector::sample(double x) // Sample the signal
{
    // Save the sample into the buffer
    if(bufsize > sample_cnt) buffer[sample_cnt] = x;
    sample_cnt++; // Update the number of samples detected
    if(prev < 0.0 && x > 0.0) edge_cnt++; // Number of positive edges
    prev = x; // Hold this sample in the prev variable to detect edges
}

// Get the normalized frequency
double FrequencyDetector::getNormalizedFrequency()
{
    // If we don't use the buffer, then we are simply using the counter.
    // Also if we don't have a valid pointer to a fftw calculator.
    if(bufsize < sample_cnt || calc_fftw_dft_1d == NULL)
        return (double)edge_cnt/sample_cnt;
    else // Use the fftw calculator to get the frequency
        return calc_fftw_dft_1d->getNormalizedFrequency(buffer);
}

/*****************************************************************************/
// Implementation of the TuneDetector class
// Constructor
TuneDetector::TuneDetector(int npart, int nturns): CppPyWrapper(NULL)
{
    calc_fftw_dft_1d = new FFT1D(nturns); // Allocate space for a FFT calculator
    nparticles = 0;
    reset(npart, nturns); // Simply call reset!
}

// Reset the counters and the data arrays
void TuneDetector::reset(int npart, int nturns)
{
    calc_fftw_dft_1d->reset(nturns); // Reset the FFT calculator

    if(npart != nparticles)
    {
        nparticles = npart;
        // Resize the frequency counter arrays
        fcntx_all.resize(nparticles);
        fcnty_all.resize(nparticles);
        fcntz_all.resize(nparticles);
        fcntx_first.resize(nparticles);
        fcnty_first.resize(nparticles);
        fcntz_first.resize(nparticles);

    }
    trigcnt_all = trigcnt_first = 0; // Initialize the trigger counters
    int i; // Loop over all particles

    for(i=0;i<nparticles;i++)
    {
        // reset all frequency counters
        fcntx_all[i].reset();
        fcnty_all[i].reset();
        fcntz_all[i].reset();
        fcntx_first[i].reset(nturns, calc_fftw_dft_1d);
        fcnty_first[i].reset(nturns, calc_fftw_dft_1d);
        fcntz_first[i].reset(nturns);
    }
}

// Tracking function
void TuneDetector::trackBunch(Bunch* bunch, double T[][5], bool isfirst)
{
    // Loop index and number of particles
    int i, npart_in_bunch = bunch->getSize();

    // If the number of particles in the bunch is more than the allocated size
    if(npart_in_bunch > nparticles)
        ORBIT_MPI_Finalize("monitor.tunedetector: trackBunch encounted more particles than allocated.");

    // Loop over all particles
    double xnow, ynow, znow;
    for (i = 0; i < npart_in_bunch; i++)
    {
        // Obtain betatron oscillation coordinates
        double dE = bunch->dE(i); 
        double x  = bunch->x(i)  - T[0][4]*dE;
        double xp = bunch->xp(i) - T[1][4]*dE;
        double y  = bunch->y(i)  - T[2][4]*dE;
        double yp = bunch->yp(i) - T[3][4]*dE;

        // Apply the Floquet transform matrix to the 4D data
        // to obtain the transverse momenta in normal mode space
        xnow = T[0][0]*x + T[0][1]*xp + T[0][2]*y + T[0][3]*yp;
        ynow = T[2][0]*x + T[2][1]*xp + T[2][2]*y + T[2][3]*yp;
        znow = bunch->z(i); // Pass through the z coordinate!

        // Update frequency counters for all triggers
        fcntx_all[i].sample(xnow);
        fcnty_all[i].sample(ynow);
        fcntz_all[i].sample(znow);

        // Update lattice entrance frequency counters
        if(isfirst)
        {
            fcntx_first[i].sample(xnow);
            fcnty_first[i].sample(ynow);
            fcntz_first[i].sample(znow);
        }
    }
    trigcnt_all++; // Update the total number of triggers
    if(isfirst)trigcnt_first++;// Update the number for lattice entrance
}

// Tune calculators
// General function to calculate tunes
void TuneDetector::gettunes(std::vector<double>& mu,
                std::vector<FrequencyDetector>& fcnt_all,
                std::vector<FrequencyDetector>& fcnt_first,
                bool usefirst)
{
    int nturns = fcnt_first[0].getBufferSize(); // Extract the number of turns
    // First resize the array of tunes
    mu.resize(nparticles);
    int i; // Loop over all particles
    for(i=0;i<nparticles;i++)
    {
        if(usefirst) // Are we instructed to use only the first node?
        {
            mu[i] = fcnt_first[i].getNormalizedFrequency();
        }
        // Do we have any lattice entrance triggers
        else if(fcnt_first[i].getSampleCount() == 0)
        {
            // If not, this is our best estimate
            mu[i] = (double)fcnt_all[i].getEdgeCount()/nturns;
        }
        else // Combine information from all tune detectors
        {
            // Calculate the normalized frequency using all triggers
            double normalized_freq_all = (double)fcnt_all[i].getEdgeCount()/nturns;
            // Calculate the normalized frequency using lattice entrance triggers
            double normalized_freq_first = fcnt_first[i].getNormalizedFrequency();

            // Use the normalized frequency over all triggers to obtain
            // the integer tune.
            double mu_int = std::floor(normalized_freq_all);
            // Get an estimate of the fractional tune.
            double mu_frac_estimate = normalized_freq_all - mu_int;
            // If the estimate is lower than 0.5, use the normalized frequency
            // at the lattice entrance directly
            if(mu_frac_estimate < 0.5)
                mu[i] = mu_int + normalized_freq_first;
            // Otherwise flip the frequency value
            else mu[i] = mu_int + 1.0 - normalized_freq_first;
        }

    }
}

/*****************************************************************************/
