// Implements the electron beam object

#ifndef ELECTRONBEAM_H
#define ELECTRONBEAM_H

//MPI Function Wrappers
#include "orbit_mpi.hh"
#include "wrap_mpi_comm.hh"

//pyORBIT utils
#include "CppPyWrapper.hh"

// Standard includes
#include <iostream>
#include <cmath>
#include <vector>

// Orbit constants
#include "OrbitConst.hh"

// More physical constants
#define kB 1.38064852e-23 // Boltzmann constant in J/K
#define epsilon_0 8.8541878128e-12 // Vacuum permittivity in C^2/(Nm^2)
#define mass_electron_MKS 9.109383702e-31 // Electron mass in kg
#define one_over_4pi_episilon_0 8.98755179e9 // (Nm^2)/C^2

// Electron cooler class
class ElectronBeam: public OrbitUtils::CppPyWrapper
{
    private:
        // Nominal beam velocity/c, Voltage jitter and voltage error
        double beta0, sigma_V, delta_V;
        // The beam current, pipe-radius, space charge depression
        double current, pipe_radius;
        // Derived parameters
        double Tpara, Tperp; // Temperatures of the beam
        double Vsc; // Space charge depression
        // Radius coordinate of each bin, density and z velocity distribution of electron beam
        std::vector<double> profile_r, profile_ne, profile_Er, profile_V;
        // Prepare a flat electron beam
        void prepareFlatBeam();

    public:
        int nbins; // Number of bins
        double beta, gamma, rmax; // Relativistic parameters and  maximum radius of beam

        // Constructor for uniform electron beam
        ElectronBeam(double beta_e, double Vjitter, double Tperp_e, double radius_e,
                double current_e, double r_pipe=0.0);
        // Destructor
        ~ElectronBeam();
        // Apply voltage error
        void setdeltaV(double dE);
        // Return beam density
        double getBeamDensity(double r);
        // Return beam properties at given radius
        void getBeamPropertiesatRadius(double x, double y, double B, double &nebeam, double *vebeam);
        // Get beam temperatures
        double getTpara();
        double getTperp();
        // Get local plasma frequency
        double getOmegap(double r);
        // Get local mean particle spacing
        double getMeanSpacing(double r);
        // Get the static field information at given radius
        void getStaticField(double r, double &Er, double &V);
};

#endif // End of header
