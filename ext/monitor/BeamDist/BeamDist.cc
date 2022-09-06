// Implements the beam distribution monitor

#include "BeamDist.hh"

using namespace OrbitUtils;

/*****************************************************************************/
// Constructor
BeamDist::BeamDist(int npart, TRIG_TYPE trig_mode): CppPyWrapper(NULL)
{
    // std::cout << "Creating an instance of BeamDist with npart = " << npart
    //     << " trig_mode = " << trig_mode << std::endl; 
    nparticles = npart; // Store the number of particles
    // Resize the particle data storage
    x.resize(nparticles); y.resize(nparticles); z.resize(nparticles);
    xp.resize(nparticles); yp.resize(nparticles); dE.resize(nparticles);
    trigger_mode = trig_mode; // Store the trigger mode
    // Initialize the trigger
    if(trigger_mode == TRIG_SINGLE) trigger_line = 0; // Not triggered by default!
    else if(trigger_mode == TRIG_RUN) trigger_line = 1; // Triggered always
}

/*****************************************************************************/
// Destructor
BeamDist::~BeamDist()
{
    return; // Nothing to do here!
}

/*****************************************************************************/
// Trigger acquisition
void BeamDist::trigger()
{
    trigger_line = 1; // Just set the trigger to high
}

/*****************************************************************************/
// Tracking function
void BeamDist::trackBunch(Bunch* bunch)
{
    if(!trigger_line)
    {
        // std::cout << "Not triggered!" << std::endl;
        return; // Acquision is not triggered. Do nothing!
    }

    // std::cout << "Triggered!" << std::endl;

    // Loop index and number of particles
    int i, npart_in_bunch = bunch->getSize();

    // If the number of particles in the bunch don't match our allocated size
    // then resize the vectors. Although this can be time consuming. Should
    // come up with a clever way to deal with beam loss.
    if(npart_in_bunch != nparticles)
    {
        nparticles = npart_in_bunch; // Set the new size
        // Resize the particle data storage
        x.resize(nparticles); y.resize(nparticles); z.resize(nparticles);
        xp.resize(nparticles); yp.resize(nparticles); dE.resize(nparticles);
    }

    // Loop over all particles
    for (i = 0; i < nparticles; i++)
    {
        // Store phase space coordinates of all particles in bunch
        x[i] = bunch->x(i); y[i] = bunch->y(i); z[i] = bunch->z(i);
        xp[i] = bunch->xp(i); yp[i] = bunch->yp(i); dE[i] = bunch->pz(i);
    }
    
    // Reset the trigger if in single shot mode
    if(trigger_mode == TRIG_SINGLE) trigger_line = 0;
}
