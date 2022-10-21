// Implements the tune detector

#include "TuneDetector.hh"

using namespace OrbitUtils;

/*****************************************************************************/
// Constructor
TuneDetector::TuneDetector(int npart): CppPyWrapper(NULL)
{
    nparticles = 0;
    reset(npart); // Simply call reset!
}

/*****************************************************************************/
// Destructor
TuneDetector::~TuneDetector()
{
    return; // Nothing to do here!
}

/*****************************************************************************/
// Reset the counters and the data arrays
void TuneDetector::reset(int npart)
{
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
        fcntx_first[i].reset();
        fcnty_first[i].reset();
        fcntz_first[i].reset();
    }
}

/*****************************************************************************/
// Tracking function
void TuneDetector::trackBunch(Bunch* bunch, double T[][4], bool isfirst)
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
        // Obtain normalized coords
        double x = bunch->x(i), y = bunch->y(i),
               xp = bunch->xp(i), yp = bunch->yp(i);
        // Apply the Floquet transform matrix to the 4D data
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

/*****************************************************************************/
// Tune calculators
// General function to calculate tunes
void TuneDetector::gettunes(std::vector<double>& mu,
                std::vector<FrequencyCounter>& fcnt_all,
                std::vector<FrequencyCounter>& fcnt_first, int nturns)
{
    // First resize the array of tunes
    mu.resize(nparticles);
    int i; // Loop over all particles
    for(i=0;i<nparticles;i++)
    {
        // Do we have any lattice entrance triggers
        if(fcnt_first[i].getSampleCount() == 0)
        {
            // If not, this is our best estimate
            mu[i] = (double)fcnt_all[i].getEdgeCount()/nturns;
        }
        else
        {
            // Calculate the normalized frequency using all triggers
            double normalized_freq_all = (double)fcnt_all[i].getEdgeCount()/
                fcnt_first[i].getSampleCount();
            // Calculate the normalized frequency using lattice entrance triggers
            double normalized_freq_first = (double)fcnt_first[i].getEdgeCount()/
                fcnt_first[i].getSampleCount();

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

