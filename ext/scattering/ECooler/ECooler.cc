// Implements electron cooling

#include "ECooler.hh"

using namespace OrbitUtils;

/*****************************************************************************/
// Constructor - for dummy cooling
ECooler::ECooler(double damp_x, double damp_y, double damp_z): CppPyWrapper(NULL)
{
    std::cout << WARNING_PREFIX << "Using a dummy cooler model." << std::endl;
    cooler_model = ECOOLER_DUMMY;
    dx = damp_x;
    dy = damp_y;
    dz = damp_z;
}

// Constructor for all internal models
ECooler::ECooler(ElectronBeam *ebeam, double lmax, ECOOLER_MODEL model, double B,
        int Z, double fscale): CppPyWrapper(NULL)
{
    max_cooler_seg_length = lmax; // Set the maximum length of a cooler segment
    Bcooler = B; // Set the magnetic field
    Zions = Z; // Charge of ions
    force_scaling = fscale; // Set the arbitrary force scaling
    CoolingBeam = ebeam; // Set the cooling beam
    cooler_model = model; // Set the cooler model
    include_cooling_force = include_electron_lens_kick = true; 
    // Do necessary preparation for different models
    switch(cooler_model)
    {
        case ECOOLER_NON_MAGNETIC:
            createNonMagneticLookupTable();
            break;
        case ECOOLER_PARKHOMCHUK:
            prepareParkhomchuk();
            break;
        default:
            ORBIT_MPI_Finalize("Cooler model not recognized.");
    }
    // Prepare default rotation and boost matrices
    prepareRotationMatrices(0.0,0.0);
    prepareBoostMatrices();

}

// Constructor for LUT model
ECooler::ECooler(ElectronBeam *ebeam, double prefactor, std::vector<double> &vz,
                 std::vector<double> &vr, std::vector<double> &fpara,
                 std::vector<double> &fperp, int Z): CppPyWrapper(NULL)
{
    cooler_model = ECOOLER_LUT; // Set the cooler model type
    lut_force_prefactor = prefactor; // Prefactor for LUT force
    force_scaling = 1.0; // Set the default scale factor to 1
    CoolingBeam = ebeam; // Set the cooling beam
    Zions = Z; // Charge of ions

    include_cooling_force = include_electron_lens_kick = true;

    if(CoolingBeam->nbins > 1)
        std::cout << WARNING_PREFIX << "Local plasma parameters ignored"
            " for non-uniform electron beam." << std::endl;
    // Copy the LUT velocity arrays from input
    lut_vz = vz; lut_vr = vr;
    int i, j, zpts = lut_vz.size(), rpts = lut_vr.size();

    lut_vz_max = *std::max_element(lut_vz.begin(), lut_vz.end());
    lut_vr_max = *std::max_element(lut_vr.begin(), lut_vr.end());

    // Allocate 2D GSL interpolation objects
    force_gsl_zacc = gsl_interp_accel_alloc();
    force_gsl_racc = gsl_interp_accel_alloc();
    fpara_gsl_spline2d = gsl_spline2d_alloc(gsl_interp2d_bicubic, zpts, rpts);
    fperp_gsl_spline2d = gsl_spline2d_alloc(gsl_interp2d_bicubic, zpts, rpts);
    lut_fpara.resize(zpts*rpts); // Resize the fpara array
    lut_fperp.resize(zpts*rpts); // Resize the fperp array

    for(i=0; i<zpts; i++) // Loop over all indices
        for(j=0; j<rpts; j++)
        {
            // Set the force values
            gsl_spline2d_set(fpara_gsl_spline2d, lut_fpara.data(), i, j, fpara[zpts*j+i]);
            gsl_spline2d_set(fperp_gsl_spline2d, lut_fperp.data(), i, j, fperp[zpts*j+i]);
        }

    // Initialize interpolation
    gsl_spline2d_init(fpara_gsl_spline2d, lut_vz.data(), lut_vr.data(),
            lut_fpara.data(), zpts, rpts);
    gsl_spline2d_init(fperp_gsl_spline2d, lut_vz.data(), lut_vr.data(),
            lut_fperp.data(), zpts, rpts);

    // Prepare default rotation and boost matrices
    prepareRotationMatrices(0.0,0.0);
    prepareBoostMatrices();
}

// Destructor
ECooler::~ECooler()
{
    // Explicitly clear GSL lookup table data if allocated
    if(fpara_gsl_spline2d)gsl_spline2d_free(fpara_gsl_spline2d);
    if(fperp_gsl_spline2d)gsl_spline2d_free(fperp_gsl_spline2d);
    if(force_gsl_zacc)gsl_interp_accel_free(force_gsl_zacc);
    if(force_gsl_racc)gsl_interp_accel_free(force_gsl_racc);
}

// Create the non-magnetic cooling lookup table
void ECooler::createNonMagneticLookupTable()
{
    int mpi_size=1, mpi_rank=0, mpi_init_flag, i, j, rpts, zpts;
    double vr, vz, d_e, sigma_v_sqr, bmin, bmax, fpara,
           fperp, send_data, recv_data, beta_e, gamma_e;
    if(CoolingBeam->nbins > 1)
        std::cout << WARNING_PREFIX << "Local plasma parameters ignored"
            " for non-uniform electron beam." << std::endl;
    // First gather important parameters of the electron beam
    beta_e = CoolingBeam->beta;
    gamma_e = pow(1-beta_e*beta_e, -0.5);
    sigma_v_e_z = sqrt(CoolingBeam->getTpara()*kB/mass_electron_MKS);
    sigma_v_e_r = sqrt(CoolingBeam->getTperp()*kB/mass_electron_MKS);
    sigma_v_sqr = sigma_v_e_r*sigma_v_e_r + sigma_v_e_z*sigma_v_e_z;
    omega_p = CoolingBeam->getOmegap(CoolingBeam->rmax/2.0);
    d_e = CoolingBeam->getMeanSpacing(CoolingBeam->rmax/2.0);
    nm_force_prefactor = force_scaling*4.0*M_PI*mass_electron_MKS*r_ec_sqr*r_ec_sqr/
        (sigma_v_e_r*sigma_v_e_r);

    bmax = std::min(std::max(sqrt(sigma_v_sqr)/omega_p, d_e),
            max_cooler_seg_length/(beta_e*gamma_e));
    bmin = r_ec_sqr/sigma_v_sqr; // Minimum impact parameter
    nm_force_asymptotic_Lc = log(bmax/bmin);
    if(sigma_v_e_z > sigma_v_e_r)
        std::cout << WARNING_PREFIX << "The longitudinal temperature of the"
            " electron beam exceeds transverse temperature." << std::endl;

    // Integration parameters
    double params[] = {0.0, 0.0, 2.0*sigma_v_e_z*sigma_v_e_z/(sigma_v_e_r*sigma_v_e_r),
        2.0, sigma_v_e_z*sigma_v_e_z/(sigma_v_e_r*sigma_v_e_r)+2.0,
        sigma_v_e_z*pow(2.*M_PI,1.5)/sigma_v_e_r,
        abs(Zions)*r_ec_sqr/(sigma_v_e_r*sigma_v_e_r), omega_p/sigma_v_e_r, d_e,
        max_cooler_seg_length/(beta_e*gamma_e), sigma_v_e_z/sigma_v_e_r};

    // Declare integration limits - won't integrate to \infty
    dblvec ll {0.0, -3.0*sigma_v_e_z/sigma_v_e_r, 0.0};
    dblvec ul {3.0, 3.0*sigma_v_e_z/sigma_v_e_r, M_PI};
    // Correction factor to approximately scale the force
    // Int = [1-\exp(-r_0^2/\sigma_r^2)]*\erf(z_0/(\sqrt{2}\sigma_z))
    double correction_factor = (1.0-exp(-4.5))*erf(3.0/sqrt(2.0));


    // Set up fpara integrator
    ndim_integrator intfpara(3, ll, ul);
    intfpara.set_integrand(0, nm_force_integrand_r);
    intfpara.set_integrand(1, nm_force_integrand_z);
    intfpara.set_integrand(2, nm_fpara_integrand_phi);
    intfpara.set_integrator_type(1, NDIM_INT_QAGP);
    intfpara.set_singularity(1, nm_force_singularity_z);
    intfpara.set_verbosity(1);
    intfpara.set_name("FPARA");

    // Set up fperp integrator
    ndim_integrator intfperp(3, ll, ul);
    intfperp.set_integrand(0, nm_force_integrand_r);
    intfperp.set_integrand(1, nm_force_integrand_z);
    intfperp.set_integrand(2, nm_fperp_integrand_phi);
    intfperp.set_integrator_type(0, NDIM_INT_QAGP);
    intfperp.set_singularity(0, nm_force_singularity_r);
    intfperp.set_integrator_type(1, NDIM_INT_QAGP);
    intfperp.set_singularity(1, nm_force_singularity_z);
    intfperp.set_verbosity(1);
    intfperp.set_name("FPERP");

    // Create sampling arrays
    create_nonuniform_sampling(lut_vz, NM_LUT_VELOCITY_MAX*sigma_v_e_r/sigma_v_e_z);
    create_nonuniform_sampling(lut_vr, NM_LUT_VELOCITY_MAX);
    zpts = lut_vz.size(); rpts = lut_vr.size();

    // Set up MPI communications
    ORBIT_MPI_Initialized(&mpi_init_flag);
    if(mpi_init_flag)
    {
        ORBIT_MPI_Comm_size(MPI_COMM_WORLD,  &mpi_size);
        ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    }

    // Start the LUT calculations
    if(!mpi_rank)std::cout << WARNING_PREFIX << "Creating lookup table for"
        " non-magnetic cooling. This will take a while." << std::endl;
    auto tic = std::chrono::high_resolution_clock::now();

    // Allocate 2D GSL interpolation objects
    force_gsl_zacc = gsl_interp_accel_alloc();
    force_gsl_racc = gsl_interp_accel_alloc();
    fpara_gsl_spline2d = gsl_spline2d_alloc(gsl_interp2d_bicubic, zpts, rpts);
    fperp_gsl_spline2d = gsl_spline2d_alloc(gsl_interp2d_bicubic, zpts, rpts);
    lut_fpara.resize(zpts*rpts); // Resize the fpara array
    lut_fperp.resize(zpts*rpts); // Resize the fperp array
   
    if(!mpi_rank) std::cout << "Calculating lut_vz = "; 
    // Calculate the cooling forces for selected points on the interpolation grid
    for(i=mpi_rank; i<zpts; i+=mpi_size) // Basically alternate between ranks
    {
        std::cout << lut_vz[i] << ", "; // Print the current z value
        for(j=0; j<rpts; j++)
        {
            // Set up integrator parameters
            params[0] = lut_vz[i]*sigma_v_e_z/sigma_v_e_r; params[1] = lut_vr[j];
            if(i==0 && j==0) // lut_vr[0] = lut_vz[0] = 0.0
                fpara = fperp = 0.0; // The force is zero for a particle at rest in the CM frame
            else // Otherwise integrate
            {
                if(i==0) fpara = 0.0;
                else fpara = intfpara.evaluate(params)/correction_factor;

                if(j==0) fperp = 0.0;
                else fperp = intfperp.evaluate(params)/correction_factor; // Needs to be debugged!!

                /*if(params[1] < 0.01)
                {
                    intfperp.set_integrand(2, nm_fperp_integrand_phi_alt); // Approximation
                    fperp = intfperp.evaluate(params)/correction_factor;
                    intfperp.set_integrand(2, nm_fperp_integrand_phi); // Reset approximation
                }
                else fperp = intfperp.evaluate(params)/correction_factor; // Needs to be debugged!!
                */
            }
            // Set the force values
            gsl_spline2d_set(fpara_gsl_spline2d, lut_fpara.data(), i, j, fpara);
            gsl_spline2d_set(fperp_gsl_spline2d, lut_fperp.data(), i, j, fperp);
        }
    }

    // Wait for all other ranks to complete
    ORBIT_MPI_Barrier(MPI_COMM_WORLD);
    if(!mpi_rank) std::cout << std::endl;

    // Get all values from other ranks! 
    for(i=0; i<zpts; i++)
        for(j=0; j<rpts; j++)
        {
            send_data = lut_fpara[i*rpts+j]; recv_data = 0.0; // Setup buffer
            // Only the value of the array elemnt in one rank will be non-zero, all others
            // will be zero. So the sum is the value of the force.
            ORBIT_MPI_Allreduce(&send_data, &recv_data, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            lut_fpara[i*rpts+j] = recv_data; // Update the array
            // Do the same with fperp
            send_data = lut_fperp[i*rpts+j]; recv_data = 0.0;
            ORBIT_MPI_Allreduce(&send_data, &recv_data, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            lut_fperp[i*rpts+j] = recv_data;
        }

    // Initialize interpolation
    gsl_spline2d_init(fpara_gsl_spline2d, lut_vz.data(), lut_vr.data(),
            lut_fpara.data(), zpts, rpts);
    gsl_spline2d_init(fperp_gsl_spline2d, lut_vz.data(), lut_vr.data(),
            lut_fperp.data(), zpts, rpts);

    // Stop measuring time and calculate the elapsed time
    auto toc = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(toc - tic);
    if(mpi_rank == 0) std::cout << "Lookup table with " << zpts*rpts <<
        " points computed in " << std::setprecision(3) << elapsed.count()*1e-6 <<
            " seconds." << std::endl;
}

// Prepare the Parkhomchuk model
void ECooler::prepareParkhomchuk()
{
    double beta_e, gamma_e, Er, Vdep, V0; 
    if(CoolingBeam->nbins > 1)
        std::cout << WARNING_PREFIX << "Local plasma parameters except density"
            " ignored for non-uniform electron beam." << std::endl;

    // Gather important parameters of the electron beam
    beta_e = CoolingBeam->beta;
    gamma_e = pow(1-beta_e*beta_e, -0.5);
    sigma_v_e_z = sqrt(CoolingBeam->getTpara()*kB/mass_electron_MKS);
    sigma_v_e_r = sqrt(CoolingBeam->getTperp()*kB/mass_electron_MKS);
    omega_p = CoolingBeam->getOmegap(CoolingBeam->rmax/2.0);
    
    // Now set up Parkhomchuk model parameters
    // $\frac{4 m_e (r_e c^2)^2}{\Delta_{e,\perp}^2}
    parkhomchuk_force_prefactor = force_scaling*4.0*mass_electron_MKS*r_ec_sqr*
        r_ec_sqr/(sigma_v_e_r*sigma_v_e_r);
    // \rho_L = \Delta_{e,\perp}/\omega_{c,e}
    parkhomchuk_rho_L = sigma_v_e_r*mass_electron_MKS/
        (OrbitConst::elementary_charge_MKS*Bcooler); // Larmor radius
    
    // Effective velocity of electrons - this will be a function of position
    // parkhomchuk_v_sqr_e_eff = sigma_v_e_z*sigma_v_e_z + Er*Er/(Bcooler*Bcooler);
    // Inverse of time of flight through cooler segment
    parkhomchuk_1_over_tau = beta_e*gamma_e*OrbitConst::c/max_cooler_seg_length;
}

// Switches for electron cooling and the e-lens kick
void ECooler::setUsageELKick(bool flag_elkick)
{
    include_electron_lens_kick = flag_elkick;
}

void ECooler::setUsageECool(bool flag_ecool)
{
    include_cooling_force = flag_ecool;
}

// Set the force scaling
void ECooler::setForceScaleFactor(double fscale)
{
    // Do necessary scaling for different models
    switch(cooler_model)
    {
        case ECOOLER_NON_MAGNETIC:
            nm_force_prefactor *= fscale/force_scaling;
            break;
        case ECOOLER_PARKHOMCHUK:
            parkhomchuk_force_prefactor *= fscale/force_scaling;
        case ECOOLER_LUT:
            lut_force_prefactor *= fscale/force_scaling;
            break;
        default:
            break;
    }
    force_scaling = fscale;
}

// Set errors for the electron beam, including offsers and misalignments
void ECooler::setErrors(double dE, double dtheta, double dphi)
{
    deltaE_e = dE; delta_theta_e = dtheta; delta_phi_e = dphi; // Offsets
    CoolingBeam->setdeltaV(deltaE_e); // Pass the energy error to the e beam
    prepareRotationMatrices(dtheta, dphi); // Make rotation matrices
    prepareBoostMatrices(); // Make boost matrices
    /*std::cout << "dE = " << dE << " eV, dx = " << dx  << " m, dy = " << dy
              << " m, dz = " << dz << " m, dtheta = " << dtheta
              << " rad, dphi = " << dphi << " rad" << std::endl;*/

    // Some debugging
    /* double vec_lab2[4], vec_pipe_e2[4], vec_e[4], vec_pipe_e[4], vec_lab[] = {1.0027,0.0,0.0,0.0732};
    PRINT4VECTOR("vec_lab", vec_lab);
    transform_vector(4, boost_PBF_from_LF, vec_lab, vec_pipe_e);
    PRINT4VECTOR("vec_pipe_e", vec_pipe_e);
    transform_vector(4, rotation_BF_from_PBF, vec_pipe_e, vec_e);
    PRINT4VECTOR("vec_e", vec_e);
    transform_vector(4, rotation_PBF_from_BF, vec_e, vec_pipe_e2);
    PRINT4VECTOR("vec_pipe_e2", vec_pipe_e2);
    transform_vector(4, boost_LF_from_PBF, vec_pipe_e2, vec_lab2);
    PRINT4VECTOR("vec_lab2", vec_lab2);*/
}

// Get errors for the electron beam, including offsers and misalignments
void ECooler::getErrors(double &dE, double &dtheta, double &dphi)
{
    dE = deltaE_e; dtheta = delta_theta_e; dphi = delta_phi_e;
}

// Prepare rotation matrices
void ECooler::prepareRotationMatrices(double dtheta, double dphi)
{
    // Prepare matrix component values
    double cos_dtheta = cos(dtheta), sin_dtheta = sin(dtheta),
           cos_dphi = cos(dphi), sin_dphi = sin(dphi);

    int i,j;

    // Rotation matrix for transforming from frame of e beam back to pipe frame
    rotation_PBF_from_BF[0][0] = 1.0;
    rotation_PBF_from_BF[1][0] = 0.0;
    rotation_PBF_from_BF[2][0] = 0.0;
    rotation_PBF_from_BF[3][0] = 0.0;

    rotation_PBF_from_BF[0][1] = 0.0;
    rotation_PBF_from_BF[1][1] = cos_dtheta*cos_dphi;
    rotation_PBF_from_BF[2][1] = cos_dtheta*sin_dphi;
    rotation_PBF_from_BF[3][1] = -sin_dtheta;

    rotation_PBF_from_BF[0][2] = 0.0;
    rotation_PBF_from_BF[1][2] = -sin_dphi;
    rotation_PBF_from_BF[2][2] = cos_dphi;
    rotation_PBF_from_BF[3][2] = 0.0;

    rotation_PBF_from_BF[0][3] = 0.0;
    rotation_PBF_from_BF[1][3] = sin_dtheta*cos_dphi;
    rotation_PBF_from_BF[2][3] = sin_dtheta*sin_dphi;
    rotation_PBF_from_BF[3][3] = cos_dtheta;

    // Rotation matrix for transforming from pipe frame to e beam frame
    // Since the rotation matrix is Hermitian, the transpose will be the
    // inverse transformation
    for(i=0;i<4;i++)
        for(j=0;j<4;j++)
            rotation_BF_from_PBF[i][j] = rotation_PBF_from_BF[j][i];
}

// Prepare boost matrices
void ECooler::prepareBoostMatrices()
{
    // Relativistic parameters
    double beta_e = CoolingBeam->beta; // Ideally oriented to the axis of pipe
    double gamma_e = pow(1-beta_e*beta_e, -0.5);

    // Veclocity components of rotated electron beam in pipe frame
    double beta_e_x = beta_e*rotation_PBF_from_BF[1][3],
           beta_e_y = beta_e*rotation_PBF_from_BF[2][3],
           beta_e_z = beta_e*rotation_PBF_from_BF[3][3];
    double beta_e_x_over_beta_e = beta_e_x/beta_e,
           beta_e_y_over_beta_e = beta_e_y/beta_e,
           beta_e_z_over_beta_e = beta_e_z/beta_e;
    
    // Boost matrix to transform from lab to e beam pipe frame
    boost_PBF_from_LF[0][0] = gamma_e;
    boost_PBF_from_LF[0][1] = -gamma_e*beta_e_x;
    boost_PBF_from_LF[0][2] = -gamma_e*beta_e_y;
    boost_PBF_from_LF[0][3] = -gamma_e*beta_e_z;

    boost_PBF_from_LF[1][0] = -gamma_e*beta_e_x;
    boost_PBF_from_LF[1][1] = 1.0+(gamma_e-1.0)*beta_e_x_over_beta_e*beta_e_x_over_beta_e;
    boost_PBF_from_LF[1][2] = (gamma_e-1.0)*beta_e_x_over_beta_e*beta_e_y_over_beta_e;
    boost_PBF_from_LF[1][3] = (gamma_e-1.0)*beta_e_x_over_beta_e*beta_e_z_over_beta_e;

    boost_PBF_from_LF[2][0] = -gamma_e*beta_e_y;
    boost_PBF_from_LF[2][1] = (gamma_e-1.0)*beta_e_y_over_beta_e*beta_e_x_over_beta_e;
    boost_PBF_from_LF[2][2] = 1.0+(gamma_e-1.0)*beta_e_y_over_beta_e*beta_e_y_over_beta_e;
    boost_PBF_from_LF[2][3] = (gamma_e-1.0)*beta_e_y_over_beta_e*beta_e_z_over_beta_e;

    boost_PBF_from_LF[3][0] = -gamma_e*beta_e_z;
    boost_PBF_from_LF[3][1] = (gamma_e-1.0)*beta_e_z_over_beta_e*beta_e_x_over_beta_e;
    boost_PBF_from_LF[3][2] = (gamma_e-1.0)*beta_e_z_over_beta_e*beta_e_y_over_beta_e;
    boost_PBF_from_LF[3][3] = 1.0+(gamma_e-1.0)*beta_e_z_over_beta_e*beta_e_z_over_beta_e;

    // Boost matrix to transform from e beam pipe frame to lab frame
    boost_LF_from_PBF[0][0] = boost_PBF_from_LF[0][0];
    boost_LF_from_PBF[0][1] = -boost_PBF_from_LF[0][1];
    boost_LF_from_PBF[0][2] = -boost_PBF_from_LF[0][2];
    boost_LF_from_PBF[0][3] = -boost_PBF_from_LF[0][3];
    boost_LF_from_PBF[1][0] = -boost_PBF_from_LF[1][0];
    boost_LF_from_PBF[1][1] = boost_PBF_from_LF[1][1];
    boost_LF_from_PBF[1][2] = boost_PBF_from_LF[1][2];
    boost_LF_from_PBF[1][3] = boost_PBF_from_LF[1][3];
    boost_LF_from_PBF[2][0] = -boost_PBF_from_LF[2][0];
    boost_LF_from_PBF[2][1] = boost_PBF_from_LF[2][1];
    boost_LF_from_PBF[2][2] = boost_PBF_from_LF[2][2];
    boost_LF_from_PBF[2][3] = boost_PBF_from_LF[2][3];
    boost_LF_from_PBF[3][0] = -boost_PBF_from_LF[3][0];
    boost_LF_from_PBF[3][1] = boost_PBF_from_LF[3][1];
    boost_LF_from_PBF[3][2] = boost_PBF_from_LF[3][2];
    boost_LF_from_PBF[3][3] = boost_PBF_from_LF[3][3];
}

// Tracking function
void ECooler::trackBunch(Bunch* bunch, double length, double xc, double yc)
{
    // Loop index and number of particles
    int i, Nparticles = bunch->getSize();

    // Ion coordinates in various frames
    double beta0_ion, gamma0_ion, rho_BF, r_LF[4], r_BF[4], p_LF_over_mc[4],
           p_PBF_over_mc[4], p_BF_over_mc[4], v_x_BF, v_y_BF, v_z_BF;

    // Forces and fields
    double fcool_x_BF, fcool_y_BF, fcool_z_BF, Esc_rho_BF, V0,
           F_BF[4], F_PBF[4], F_LF[4];

    // Ion kick parameters
    double ion_charge, kick_prefactor_xy, kick_prefactor_s;

    SyncPart* syncPart = bunch->getSyncPart(); // Get the synchronous particle

    // Gather important parameters of the electron beam and find useful values
    if(cooler_model != ECOOLER_DUMMY)
    {
        beta0_ion = syncPart->getBeta();
        gamma0_ion = syncPart->getGamma();
        ion_charge = Zions*OrbitConst::elementary_charge_MKS;
        kick_prefactor_xy = length*1e-9/(gamma0_ion*beta0_ion*beta0_ion*bunch->getMass()*
            OrbitConst::elementary_charge_MKS);
        kick_prefactor_s = syncPart->getMomentum() * beta0_ion *kick_prefactor_xy;
    }

    // Loop over all particles
    for (i = 0; i < Nparticles; i++)
    {
        if(!bunch->flag(i))continue; // Skip dead particles
        if(cooler_model == ECOOLER_DUMMY)
        {
            bunch->px(i) *= (1.0 - dx*length);
            bunch->py(i) *= (1.0 - dy*length);
            bunch->pz(i) *= (1.0 - dz*length);
        }
        else if(include_cooling_force || include_electron_lens_kick)
        {
            // Step 1: Obtain 6D ion coords in lab frame (LF).
            double x_ion = bunch->x(i), y_ion = bunch->y(i), z_ion = bunch->z(i),
                   px_ion = bunch->px(i), py_ion = bunch->py(i), pz_ion = bunch->pz(i);
            // Declare 4-momentum components of ion
            // First the transverse, p_x=p_0 x', where x' = px_ion.
            p_LF_over_mc[1] = gamma0_ion*beta0_ion*px_ion;
            p_LF_over_mc[2] = gamma0_ion*beta0_ion*py_ion;
            // Then the longitudinal, p_z=p_0+\frac{\delta E}{\beta_0 mc^2}
            // where \delta E = pz_ion
            p_LF_over_mc[3] = gamma0_ion*beta0_ion
                                + pz_ion/(beta0_ion*bunch->getMass());
            // Finally get the energy, p^0=\sqrt{p_x^2+p_y^2+p_z^2+m^2 c^2}
            p_LF_over_mc[0] = sqrt(p_LF_over_mc[1]*p_LF_over_mc[1]
                                  +p_LF_over_mc[2]*p_LF_over_mc[2]
                                  +p_LF_over_mc[3]*p_LF_over_mc[3]+1.0);


            // Step 2: Convert LF coords to beam frame (BF) coords
            // Step 2a: Calculate ion transverse position.
            // We always start at z=0,t=0 inside each element.
            r_LF[0] = r_LF[3] = 0.0; r_LF[1] = x_ion-xc; r_LF[2] = y_ion-yc;
            // Rotate LF coords into the electron beam orientation.
            // Note that we don't apply a boost since the transverse size of
            // the electron beam is measured in the lab frame. Also note
            // that if the ion beam and the electron beams are well aligned
            // then the transverse size of the beam does not change between
            // the lab frame and the beam frame.
            transform_vector(4, rotation_BF_from_PBF, r_LF, r_BF);
            rho_BF = sqrt(r_BF[1]*r_BF[1]+r_BF[2]*r_BF[2]); // The radial pos
            // Step 2b: Calculate 4-momentum in BF
            // First boost the LF coordinates to the 'pipe beam frame' (PBF).
            // PBF is aligned to the lab frame but travels at the same
            // velocity as the electron beam.
            transform_vector(4, boost_PBF_from_LF, p_LF_over_mc, p_PBF_over_mc);
            // Then rotate the boosted 4-momentum to align with the e beam.
            transform_vector(4, rotation_BF_from_PBF, p_PBF_over_mc, p_BF_over_mc);
            // Step 2c: Calulate velocities in BF
            v_x_BF = OrbitConst::c*p_BF_over_mc[1]/p_BF_over_mc[0];
            v_y_BF = OrbitConst::c*p_BF_over_mc[2]/p_BF_over_mc[0];
            v_z_BF = OrbitConst::c*p_BF_over_mc[3]/p_BF_over_mc[0];
           

            // Step 3: Calculate the net force on a single ion
            // Set all forces to 0 by default
            fcool_x_BF = fcool_y_BF = fcool_z_BF = Esc_rho_BF = 0.0;
            // Get the cooling force
            if(include_cooling_force)
                getcoolingforce(r_BF[1], r_BF[2], v_x_BF, v_y_BF, v_z_BF,
                        fcool_x_BF, fcool_y_BF, fcool_z_BF);
            // Get the radial electric field due to the space charge of the e beam
            if(include_electron_lens_kick)
                CoolingBeam->getStaticField(rho_BF, Esc_rho_BF, V0);
            
            // Calculate the 4-force on the ion in BF.
            // F^\mu \equiv [\gamma \vec{f} \cdot \vec{\beta}, \gamma \vec{f}]
            F_BF[1] = p_BF_over_mc[0]*(fcool_x_BF + ion_charge*Esc_rho_BF*r_BF[1]/rho_BF);
            F_BF[2] = p_BF_over_mc[0]*(fcool_y_BF + ion_charge*Esc_rho_BF*r_BF[2]/rho_BF);
            F_BF[3] = p_BF_over_mc[0]*fcool_z_BF;
            F_BF[0] = (F_BF[1]*v_x_BF + F_BF[2]*v_y_BF + F_BF[3]*v_z_BF)/OrbitConst::c;
            // Now transform the 4-force to the lab frame
            transform_vector(4, rotation_PBF_from_BF, F_BF, F_PBF);
            transform_vector(4, boost_LF_from_PBF, F_PBF, F_LF);


            // Step 4: Apply kicks!
            bunch->px(i) += kick_prefactor_xy*F_LF[1]/p_LF_over_mc[0];
            bunch->py(i) += kick_prefactor_xy*F_LF[2]/p_LF_over_mc[0];
            bunch->pz(i) += kick_prefactor_s* F_LF[3]/p_LF_over_mc[0];
        }
    }
}

// Get the cooling forces
int ECooler::getcoolingforce(double x, double y, double vx, double vy, double vz,
        double &fx, double &fy, double &fz)
{
    int success;
    double prefactor, n_e, v_e[3];
    // Get the electron beam density and mean relative velocity in the beam frame
    CoolingBeam->getBeamPropertiesatRadius(x, y, Bcooler, n_e, v_e);
    // To make sure we can model dragging of ion beam energy we convert the
    // longitudinal ion velocity to the rest frame of the electron beam at the
    // radius it's interacting with.
    vz -= v_e[2];
    double r = sqrt(x*x+y*y), visqr, vr = sqrt(vx*vx+vy*vy), abs_vz = abs(vz),
           fsign = (vz < 0.0) - (vz > 0.0);
    double fpara = GSL_NAN, fperp = GSL_NAN; // For non-magnetic model and LUT model
    double parkhomchuk_v_sqr_e_eff, bmax, bmin, Lp, vsqr, vsqr_3_2; // For Parkhomchuk model
    fx = fy = fz = 0.0; // By default force is 0 
    if(r > CoolingBeam->rmax)
        return 1; // We're outside the beam, cooling force is 0
    switch(cooler_model)
    {
        case ECOOLER_NON_MAGNETIC:
            prefactor = Zions*Zions*nm_force_prefactor * n_e;
            // Normalized ion velocity square
            visqr = (vz*vz + vr*vr)/(sigma_v_e_r*sigma_v_e_r);
            // Use the lookup table within certain bounds
            if(visqr < NM_LUT_VELOCITY_MAX*NM_LUT_VELOCITY_MAX)
            {
                fpara = fsign*prefactor*gsl_spline2d_eval(fpara_gsl_spline2d,
                        abs_vz/sigma_v_e_z, vr/sigma_v_e_r, force_gsl_zacc,
                        force_gsl_racc);
                fperp = -prefactor*gsl_spline2d_eval(fperp_gsl_spline2d,
                        abs_vz/sigma_v_e_z, vr/sigma_v_e_r, force_gsl_zacc,
                        force_gsl_racc);
            }
            // If any of the forces are nan, then use an asymptote
            if(isnan(fpara))
               fpara =  -prefactor*nm_force_asymptotic_Lc*
                   (vz/sigma_v_e_r)/pow(visqr, 1.5);
            if(isnan(fperp))
               fperp =  -prefactor*nm_force_asymptotic_Lc*
                   (vr/sigma_v_e_r)/pow(visqr, 1.5);
            success = 1; // We were successful
            break;

        case ECOOLER_PARKHOMCHUK:
            prefactor = Zions*Zions*parkhomchuk_force_prefactor * n_e 
                *sigma_v_e_r*sigma_v_e_r;
            visqr = vz*vz + vr*vr; // Ion velocity square
            // Effective velocity of electron beam
            parkhomchuk_v_sqr_e_eff = v_e[0]*v_e[0]+v_e[1]*v_e[1]+sigma_v_e_z*sigma_v_e_z;
            vsqr = visqr+parkhomchuk_v_sqr_e_eff;
            vsqr_3_2 = pow(vsqr, 1.5);
            bmin = Zions*r_ec_sqr/visqr; // Minimum impact parameter
            // Maximum impact parameter
            bmax = sqrt(visqr)/(parkhomchuk_1_over_tau+omega_p);
            // Coulomb log
            Lp = log((bmax+bmin+parkhomchuk_rho_L)/(bmin+parkhomchuk_rho_L));
            // Cooling forces
            fpara = -prefactor*vz*Lp/vsqr_3_2;
            fperp = -prefactor*vr*Lp/vsqr_3_2;
            success = 1; // We were successful
            break;

        case ECOOLER_LUT:
            prefactor = Zions*Zions*lut_force_prefactor * n_e;
            // Use the lookup table within the bounds of the domain
            if(abs_vz < lut_vz_max && vr < lut_vr_max)
            {
                fpara = fsign*prefactor*gsl_spline2d_eval(fpara_gsl_spline2d,
                            abs_vz, vr, force_gsl_zacc, force_gsl_racc);
                fperp = -prefactor*gsl_spline2d_eval(fperp_gsl_spline2d,
                            abs_vz, vr, force_gsl_zacc, force_gsl_racc);
            }
            if(isnan(fpara))fpara=0.0; // Out of the domain values are set to 0
            if(isnan(fperp))fperp=0.0;
            success = 1; // We were successful
            break;

        default:
            success = 0;
    }
    if(success)
    {
        fx = fperp*vx/vr;
        fy = fperp*vy/vr;
        fz = fpara;
    }
    return success; // Return success flag
}

// Copy the Look-Up Table
int ECooler::copyLUT(double &prefactor, std::vector<double> &vz, std::vector<double> &vr,
                     std::vector<double> &fpara, std::vector<double> &fperp)
{
    // Only two cooler models use the LUT
    if(cooler_model != ECOOLER_NON_MAGNETIC && cooler_model != ECOOLER_LUT)
        return 0; // Nothing to copy
    // Copy the force
    fpara = lut_fpara; fperp = lut_fperp;
    if(cooler_model == ECOOLER_NON_MAGNETIC)
    {
        prefactor = nm_force_prefactor;
        for(double val: lut_vz) vz.push_back(sigma_v_e_z*val);
        for(double val: lut_vr) vr.push_back(sigma_v_e_r*val);
    }
    if(cooler_model == ECOOLER_LUT)
    {
        prefactor = lut_force_prefactor;
        vz = lut_vz; vr = lut_vr;
    }
    return 1;
}

/*****************************************************************************/
// Support functions
// Electron cooling integration
// xvec = {\phi, v_{e,z}, v_{e,r}}
// params = {v_{i,z}, v_{i,r}, 2\sigma_{v_{e,z}}^2, 2\sigma_{v_{e,r}}^2,
// \sigma_{v_{e,z}}^2 + 2\sigma_{v_{e,r}}^2, (2\pi)^{3/2} \sigma_{v_{e,z}},
// 2Z r_e c^2, \omega_p, d_e, l_\text{cooler}, \sigma_{v_{e,z}}}
// params and xvec are normalized with \sigma_{v_{e,r}}

// Non-magnetic longitudinal cooling force phi integrand
// Integrated from 0 to \pi
double nm_fpara_integrand_phi(double *xvec, void *params_void)
{
    double *params = (double *) params_void; // Import the params as double
    // Get the integration variables
    double phi = xvec[2], vez = xvec[1], ver = xvec[0];
    double diff_z = params[0] - vez; // This should help numerically
    if(abs(diff_z)<1e-30) return 0.0; // Return 0!
    // Square of relative velocity
    double vrelsqr = diff_z*diff_z + 2.0*params[1]*ver*cos(phi) +
        ver*ver + params[1]*params[1];
    // Range of impact parameters
    double rhomin = params[6]/vrelsqr;
    double rhomax = std::min(std::max(sqrt(vrelsqr+params[4])/params[7],
                params[8]),params[9]);
    double Lc = log(rhomax/rhomin); // Coulomb log
    if(Lc < 0.0) return 0.0;
    else return 2.0*Lc*diff_z/pow(vrelsqr, 1.5);
}

// Non-magnetic transverse cooling force phi integrand
// Integrated from 0 to \pi
double nm_fperp_integrand_phi(double *xvec, void *params_void)
{
    double *params = (double *) params_void; // Import the params as double
    // Get the integration variables
    double phi = xvec[2], vez = xvec[1], ver = xvec[0];
    double diff_r = params[1] - ver*cos(phi);
    double diff_z = params[0] - vez;
    if(abs(diff_r)<1e-30) return 0.0; // Return 0!
    // Square of relative velocity
    double vrelsqr = diff_z*diff_z - 2.0*params[1]*ver*cos(phi) +
        ver*ver + params[1]*params[1];
    // Range of impact parameters
    double rhomin = params[6]/vrelsqr;
    double rhomax = std::min(std::max(sqrt(vrelsqr+params[4])/params[7],
                params[8]),params[9]);
    double Lc = log(rhomax/rhomin); // Coulomb log
    if(Lc < 0.0) return 0.0;
    else return 2.0*Lc*diff_r/pow(vrelsqr, 1.5);
}

// Non-magnetic transverse cooling force phi alternate integrand
// Integrated from 0 to \pi
double nm_fperp_integrand_phi_alt(double *xvec, void *params_void)
{
    double *params = (double *) params_void; // Import the params as double
    // Get the integration variables
    double phi = xvec[2], vez = xvec[1], ver = xvec[0];
    double diff_z = params[0] - vez;
    // Square of relative velocity
    double vrelsqr = diff_z*diff_z - 2.0*params[1]*ver*cos(phi) +
        ver*ver + params[1]*params[1];
    if(vrelsqr < 1e-30) return GSL_POSINF;
    // Range of impact parameters
    double rhomin = params[6]/vrelsqr;
    double rhomax = std::min(std::max(sqrt(vrelsqr+params[4])/params[7],
                params[8]),params[9]);
    double Lc = log(rhomax/rhomin); // Coulomb log
    if(Lc < 0.0) return 0.0;
    else return 2.0*Lc*params[1]/pow(vrelsqr, 1.5);
}

// Non-magnetic cooling force z integrand
// Integrated from -3\sigma_{v_{e,z}} to 3\sigma_{v_{e,z}}
double nm_force_integrand_z(double *xvec, void *params_void)
{
    double *params = (double *) params_void; // Import the params as double
    double vez = xvec[1]; // Get the integration variables
    return exp(-vez*vez/params[2])/params[5];
}

// Non-magnetic cooling force r integrand
// Integrated from 0 to 3\sigma_{v_{e,r}}
double nm_force_integrand_r(double *xvec, void *params_void)
{
    double *params = (double *) params_void; // Import the params as double
    double ver = xvec[0]; // Get the integration variables
    return ver*exp(-ver*ver/params[3]);
}

// Non-magnetic cooling force z singularity
// Singularity at v_{i,z}
int nm_force_singularity_z(double *xvec, void *params_void, double *spts)
{
    double *params = (double *) params_void; // Import the params as double
    double abs_viz = abs(params[0]);
    // Singularity only used if v_{i,z} < \sigma_{v_{e,z}} for convergence
    if(abs_viz > 1e-20 && abs_viz < 2.0*params[10])
    {
        spts[0] = params[0];
        return 1;
    }
    else return 0; // Ignore singularity otherwise
}

// Non-magnetic cooling force r singularity
// Singularity at v_{i,r}
int nm_force_singularity_r(double *xvec, void *params_void, double *spts)
{
    double *params = (double *) params_void; // Import the params as double
    double abs_vir = abs(params[1]);
    // Singularity only used if v_{i,r} < \sigma_{v_{e,r}} for convergence
    if(abs_vir > 1e-20 && abs_vir < 3.0)
    {
        spts[0] = params[1];
        return 1;
    }
    else return 0; // Ignore singularity otherwise
}

// Interpolation helpers
void create_nonuniform_sampling(std::vector<double> &xvec, double xmax)
{
    xvec.push_back(0.0);
    double x = NM_LUT_VELOCITY_LINEAR;
    while(x<xmax)
    {
        xvec.push_back(x); //0.3
        // Don't even think about touching this!
        x += 0.5*sqrt(abs(x-1.0))+0.01;
    }
}

// Simple coordinate transformation
void transform_vector(int n, double T[][4], double *x, double *y)
{
    int i,j;
    for(i=0;i<n;i++)
    {
        y[i] = 0.0;
        for(j=0;j<n;j++) y[i] += T[i][j]*x[j];
    }
}
