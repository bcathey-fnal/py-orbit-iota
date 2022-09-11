// Implements electron cooling

#ifndef ECOOLER_H
#define ECOOLER_H

// Electron beam used to cool
#include "ElectronBeam.hh"

// Multi-dimensional integral
#include "ndim_integrator.hh"

// Interpolation from GSL
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

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
#include <iomanip>
#include <vector>
#include <chrono>
#include <algorithm>

// Useful macros
#define WARNING_PREFIX "Warning from line " << __LINE__ << " in file " << __FILE__ << ": "
#define PRINT4VECTOR(NAME,VEC) std::cout << NAME << "[] = {" << VEC[0] << ", " << VEC[1] << ", " << VEC[2] << ", " << VEC[3] << "}" << std::endl


// Hard coded parameters
// Physical constants
#define kB 1.38064852e-23 // Boltzmann constant in J/K
#define r_ec_sqr (OrbitConst::classicalRadius_electron*OrbitConst::c*OrbitConst::c)
// Non magnetic cooling look up table parameters, normalized to \sigma_v
#define NM_LUT_VELOCITY_LINEAR 0.001 // Normalized velocity used for linear force estimates
#define NM_LUT_VELOCITY_MAX 6.00 // Maximum normalized velocity for interpolation

// Define some constants to decide cooler type
typedef enum {ECOOLER_DUMMY, ECOOLER_NON_MAGNETIC, ECOOLER_PARKHOMCHUK, ECOOLER_LUT} ECOOLER_MODEL;

/*****************************************************************************/
// Electron cooler class
class ECooler: public OrbitUtils::CppPyWrapper
{
    private:
        ECOOLER_MODEL cooler_model; // Which model to use
        // Some general parameters
        double dx, dy, dz, max_cooler_seg_length, Bcooler, Zions, force_scaling;
        // Electron beam parameters
        double sigma_v_e_z, sigma_v_e_r, omega_p;
        // Offsets of the electron beam
        double deltaE_e, delta_theta_e, delta_phi_e;
        // Flags to turn on or off the cooling and the electron lens force
        bool include_electron_lens_kick, include_cooling_force;

        // Non-magnetic model parameters
        double nm_force_prefactor, nm_force_asymptotic_Lc;
        // Parkhomchuk model parameters
        double parkhomchuk_force_prefactor, parkhomchuk_v_sqr_e_eff,
               parkhomchuk_rho_L, parkhomchuk_1_over_tau;

        // LUT model data
        double lut_force_prefactor, lut_vz_max, lut_vr_max;

        // All LUT data - both external and also non-magnetic
        std::vector<double> lut_fpara;
        std::vector<double> lut_fperp;
        std::vector<double> lut_vz;
        std::vector<double> lut_vr;
        gsl_interp_accel *force_gsl_zacc = nullptr;
        gsl_interp_accel *force_gsl_racc = nullptr;
        gsl_spline2d *fpara_gsl_spline2d = nullptr;
        gsl_spline2d *fperp_gsl_spline2d = nullptr;

        // Rotation and Boost matrices
        double rotation_BF_from_PBF[4][4], rotation_PBF_from_BF[4][4],
            boost_PBF_from_LF[4][4], boost_LF_from_PBF[4][4];


        // Prepare the non-magnetic cooler model
        void createNonMagneticLookupTable();
        // Prepare the Parkhomchuk model
        void prepareParkhomchuk();

        // Prepare rotation matrices
        void prepareRotationMatrices(double dtheta, double dphi);
        // Prepare boost matrices
        void prepareBoostMatrices();

    protected:
        ElectronBeam *CoolingBeam;

    public:
        // Constructor - for dummy model
        ECooler(double damp_x, double damp_y, double damp_z);
        // Contructor for NON_MAGNETIC and PARKHOMCHUK
        ECooler(ElectronBeam *ebeam, double lmax, ECOOLER_MODEL alg, double B=0.0,
                int Z=1, double fscale=1.0);
        // Constructor for LUT
        ECooler(ElectronBeam *ebeam, double prefactor, std::vector<double> &vz,
                std::vector<double> &vr, std::vector<double> &fpara,
                std::vector<double> &fperp, int Z=1);
        // Destructor
        ~ECooler();
        // Set e-lens kick switch
        void setUsageELKick(bool flag_elkick=true);
        // Set electron cooling switch
        void setUsageECool(bool flag_ecool=true);
        // Set the force scaling
        void setForceScaleFactor(double fscale=1.0);
        // Set errors for the electron beam, including offsers and misalignments
        void setErrors(double dE, double dtheta, double dphi);
        // Get errors
        void getErrors(double &dE, double &dtheta, double &dphi);
        // Copy LUT information
        int copyLUT(double &prefactor, std::vector<double> &vz, std::vector<double> &vr,
                std::vector<double> &fpara, std::vector<double> &fperp);
        // Get cooling forces
        int getcoolingforce(double r, double vz, double vr, double &fpara, double &fperp);
        // Track the bunch passing through the cooler
        void trackBunch(Bunch* bunch, double length, double xc=0.0, double yc=0.0);
};

/*****************************************************************************/
// Support functions
// Non-magnetic cooling force integration functions
double nm_fpara_integrand_phi(double *xvec, void *params_void);
double nm_fperp_integrand_phi(double *xvec, void *params_void);
double nm_fperp_integrand_phi_alt(double *xvec, void *params_void);
double nm_force_integrand_z(double *xvec, void *params_void);
double nm_force_integrand_r(double *xvec, void *params_void);
int nm_force_singularity_z(double *xvec, void *params_void, double *spts);
int nm_force_singularity_r(double *xvec, void *params_void, double *spts);

// Interpolation helpers
void create_nonuniform_sampling(std::vector<double> &xvec, double xmax);

// Simple coordinate transformation
void transform_vector(int n, double T[][4], double *x, double *y);
#endif // End of header
