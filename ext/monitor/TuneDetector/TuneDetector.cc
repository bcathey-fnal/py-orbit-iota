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
        // Resize the particle data storage
        xn.resize(nparticles); yn.resize(nparticles); z.resize(nparticles);
        cntx.resize(nparticles); cnty.resize(nparticles); cntz.resize(nparticles);
    }
    trigcnt = 0; // Initialize the trigger count
    int i; // Loop over all particles and initialize count and saved values
    for(i=0;i<nparticles;i++)
    {
        cntx[i] = cnty[i] = cntz[i] = 0;
        xn[i] = yn[i] = z[i] = std::numeric_limits<double>::max();
    }
}

/*****************************************************************************/
// Tracking function
void TuneDetector::trackBunch(Bunch* bunch, double T[][4])
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
        double x = bunch->x(i), y = bunch->y(i), xp = bunch->xp(i), yp = bunch->yp(i);
        // Apply the Floquet transform matrix to the 4D data
        xnow = T[0][0]*x + T[0][1]*xp + T[0][2]*y + T[0][3]*yp;
        ynow = T[2][0]*x + T[2][1]*xp + T[2][2]*y + T[2][3]*yp;
        znow = bunch->z(i); // Pass through the z coordinate!
        //if(i==0)
        //    std::cout << "xn[0] = " << xn[0] << ", xnow = " << xnow << ", cntx[0] = " << cntx[0] << std::endl;
        //    std::cout << xnow << " " << T[1][0]*x + T[1][1]*xp + T[1][2]*y + T[1][3]*yp << std::endl;
        // Detect positive edge zero crossings only and update saved values
        if(xn[i] < 0.0 && xnow > 0.0) cntx[i]++; xn[i] = xnow;
        if(yn[i] < 0.0 && ynow > 0.0) cnty[i]++; yn[i] = ynow;
        if(z[i] < 0.0  && znow > 0.0) cntz[i]++; z[i]  = znow;
        
    }
    trigcnt++; // Update trigger count
}

