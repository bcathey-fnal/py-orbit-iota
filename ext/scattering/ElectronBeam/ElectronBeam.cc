// Implements electron beam for cooling

#include "ElectronBeam.hh"

using namespace OrbitUtils;

// Constructor
ElectronBeam::ElectronBeam(double beta_e, double Vjitter, double Tperp_e,
        double radius_e, double current_e, double r_pipe)
{
    // First copy the basic properties
    beta0 = beta = beta_e;
    sigma_V = Vjitter;
    Tperp = Tperp_e;
    current = current_e;
    // Set pipe radius to the beam radius if nothing is provided
    if(r_pipe > radius_e)pipe_radius = r_pipe;
    else pipe_radius = radius_e;
    rmax = radius_e; // Maximum radius of beam

    // Calculated parameters
    gamma = pow(1.0-beta*beta,-0.5); // Relativistic gamma

    // Prepare the requested beam distribution
    // TODO: Support for other shapes
    prepareFlatBeam();

}

// Destructor
ElectronBeam::~ElectronBeam()
{
    // The vectors will get deleted automatically
}

// Prepare a flat electron beam
void ElectronBeam::prepareFlatBeam()
{
    double beta_gamma = beta*gamma;
    
    // Set up the beam distribution
    nbins = 1; // Only 1 bin for uniform distribution
    // Allocate memory for distribution
    profile_r.resize(nbins); profile_ne.resize(nbins);

    // Fill in the distribution
    profile_r[0] = rmax; // The bin upper limit is just the radius of the beam
    // Density of the electron beam in the beam frame
    profile_ne[0] = current/(OrbitConst::elementary_charge_MKS*OrbitConst::c*\
            beta_gamma*M_PI*rmax*rmax);
    // Total space charge voltage depression in the beam frame
    Vsc = one_over_4pi_episilon_0*current/(OrbitConst::c*beta*gamma);
    // The longitudinal temperature in the beam frame from Voltage jitter
    Tpara = sigma_V*sigma_V*OrbitConst::elementary_charge_MKS/\
            (kB*OrbitConst::mass_electron*1e9*beta_gamma*beta_gamma);
}

// Apply voltage error
void ElectronBeam::setdeltaV(double dE)
{
    gamma = pow(1-beta0*beta0, -0.5); // Calculate reference gamma
    // Apply a change in energy
    gamma += dE*1e-9/OrbitConst::mass_electron; // Apply change in energy
    if(gamma < 1.0)
    {
        ORBIT_MPI_Finalize("ElectronBeam::setdeltaV(dE). Energy error exceeds beam voltage!");
        return;
    }
    beta = sqrt(1.0 - 1.0/(gamma*gamma));

    // Recalculate distribution data
    // TODO: what happens with other distributions?
    prepareFlatBeam();
}

// Return beam number density
double ElectronBeam::getBeamDensity(double r)
{
    if(r > rmax) return 0.0;
    else return profile_ne[0];
}

// Get electron beam properties
void ElectronBeam::getBeamPropertiesatRadius(double x, double y, double B, double &nebeam, double *vebeam)
{
    double V_at_r, Er_at_r, gamma_at_r, beta_at_r;
    double r = sqrt(x*x+y*y); // Get the radius
    vebeam[0] = vebeam[1] = vebeam[2] = 0.0;

    // If radius is larger than the upper limit of the last bin
    if(r > rmax)
    {
        // There is no beam at this radius!
        nebeam = 0.0;
        return;
    }
    if(nbins == 1) // Uniform beam, so not much to do.
    {
        getStaticField(r, Er_at_r, V_at_r); // Get the beam potential
        // Relativistic parameters of the electron beam at r
        gamma_at_r = gamma + V_at_r*1e-9/OrbitConst::mass_electron;
        beta_at_r = sqrt(1.0 - 1.0/(gamma_at_r*gamma_at_r));
        nebeam = profile_ne[0]; // Return the density
        // The longitudinal velocity transformed to the beam frame
        vebeam[2] = OrbitConst::c*(beta_at_r - beta)/(1.0 - beta_at_r*beta);
        // For non-zero radii, calculate the azimuthal velocity
        if(r > OrbitConst::tiny)
        {
            vebeam[0] = -y*Er_at_r/(r*B);
            vebeam[1] = x*Er_at_r/(r*B);
        }
        return;
    }
}

// Get temperatures
double ElectronBeam::getTpara()
{
    return Tpara;
}

double ElectronBeam::getTperp()
{
    return Tperp;
}

// Get local plasma frequency
double ElectronBeam::getOmegap(double r)
{
    double ne = getBeamDensity(r);
    return OrbitConst::elementary_charge_MKS*sqrt(ne/(mass_electron_MKS*epsilon_0));
}

// Get local mean particle spacing
double ElectronBeam::getMeanSpacing(double r)
{
    double ne = getBeamDensity(r);
    return pow(ne, -1./3.);
}

// Get the static field information at given radius
// TODO: Add Bphi due to electrons moving in the beam frame
void ElectronBeam::getStaticField(double r, double &Er, double &V)
{
    Er = V = 0.0; // By default field and potential is 0
    if(r > pipe_radius) return; // Field and potential is 0 outside pipe
    if(nbins == 1) // Uniform beam, so not much to do.
    {
        if(r > OrbitConst::tiny) // at the center of the beam the potential is 0.0
        {
            double a = rmax; // The radius of the flat beam
            if(r < a)  // If the point is within the beam
            {
                Er = -2.0*Vsc*r/(a*a); // Set the field
                V = Vsc*r*r/(a*a); // Set the potential, 0 at r=0
            }
            else
            {
                Er = -2.0*Vsc/r; // Set the field, continuous at r=a
                V = Vsc*(1.0+2.0*log(r/a)); // Set the potential, Vsc at r=a
            }
        } 
    }
}
