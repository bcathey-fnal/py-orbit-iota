// Implements the thick McMillan lens

#include "McMillan.hh"
//#define VERBOSE

using namespace OrbitUtils;

/*****************************************************************************/
// Functions for the McMillan class

// Constructor for the McMillan lens
McMillan::McMillan(double K, double A, int O, int N): CppPyWrapper(NULL)
{
    // Save parameters
    strength_over_length = K;
    width = A;
    Nslices = N;

    // Prepare symplectic integrator coefficients. Keep the order 2 for now.
    Ycoeffs = new YoshidaCoefficients(O);

    // Debug printing
#ifdef VERBOSE
    std::cout << WARNING_PREFIX << "FLoat size: " << sizeof(McDouble) << std::endl;
    std::cout << WARNING_PREFIX << "Initialized McMillan lens with K = " <<
        strength_over_length << " , A = " << width << " , N = " << Nslices << std::endl;
#endif
}

// Destructor for the McMillan lens
McMillan::~McMillan()
{
    delete Ycoeffs;
}

// Kick bunches in the McMillan lens
void McMillan::kickBunch(McDouble coeff, McDouble slice_length)
{
    int i;
    // First calculate amplitude of all kicks
    McDouble kick_prefactor = coeff*strength_over_length*slice_length;
    McDouble widthsqr = width*width;
    // Loop over all particles
    for (i = 0; i < Nparticles; i++)
    {
        if(!flag[i])continue; // Skip dead particles
        // Get particle parameters
        // x, y, z -> Position of particles (m)
        // px, py -> Momentum coordinates x' and y'
        // pz -> \delta E (GeV)
        McDouble x_i = x[i], y_i = y[i];

        McDouble r2i = x_i * x_i + y_i * y_i; // Radius coord
        McDouble delta = kick_prefactor/(r2i/widthsqr + 1.0); // McMillan kick

        McDouble dxp = delta * x_i; // x kick
        McDouble dyp = delta * y_i; // y kick

        // Finally apply changes
        px[i] += dxp;
        py[i] += dyp;
    }
}

// Transports bunch through drift (stolen from pyORBIT teapot class.)
// It can handle the negative lengths necessary for symplectic integration.
void McMillan::driftBunch(Bunch *bunch, McDouble length)
{
    McDouble KNL, phifac, dp_p; // Various parameters

    // Update the time of the synchornous particle
    SyncPart* syncPart = bunch->getSyncPart();
    McDouble v = OrbitConst::c * syncPart->getBeta();
    syncPart->setTime(syncPart->getTime() + length / v);

    // Energy parameters
    McDouble gamma2i = 1.0 / (syncPart->getGamma() * syncPart->getGamma());
    McDouble dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());

    //coordinate array [part. index][x,xp,y,yp,z,dE]

    for(int i = 0; i < Nparticles; i++)
    {
        dp_p = pz[i] * dp_p_coeff;
        KNL  = 1.0 / (1.0 + dp_p);
        x[i] += KNL * length * px[i];
        y[i] += KNL * length * py[i];
        phifac = (px[i] * px[i] + py[i] * py[i] +
                  dp_p * dp_p * gamma2i) / 2.0;
        phifac = (phifac * KNL - dp_p * gamma2i) * KNL;
        z[i] -= length * phifac;
    }
}

// Transport particles in the presence of solenoid field
void McMillan::solenoidBunch(Bunch* bunch, McDouble length, McDouble ks)
{
    // McDouble charge = bunch->getCharge();
    McDouble Bc = ks; // The effective strength
    McDouble KNL, phase, cs, sn, ck, sk, u_init, pu_init, phifac; // Various parameters

    // Update the time of the synchronous particle
    SyncPart* syncPart = bunch->getSyncPart();
    McDouble v = OrbitConst::c * syncPart->getBeta();
    syncPart->setTime( syncPart->getTime() + length/v);

    // Energy parameters
    McDouble gamma2i = 1.0 / (syncPart->getGamma() * syncPart->getGamma());
    McDouble dp_p, dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    //

    for(int i = 0; i < Nparticles; i++)
    {
        dp_p = pz[i] * dp_p_coeff;
        KNL  = 1.0 / (1.0 + dp_p);

        // No energy dependence?
        phase = Bc * length;
        //phase = KNL * Bc * length;

        // Precalculate some matrix parameters
        cs = cos(phase);
        sn = sin(phase);
        sk = sn / Bc;
        ck = (cs - 1.0) / Bc;
        u_init  =  y[i] / 2.0 + px[i]/ Bc;
        pu_init = -x[i] * Bc / 2. + py[i];

        // Perform matrix math
        McDouble xf = x[i] + sk * px[i] + ck * py[i];
        McDouble pxf = cs * px[i] - sn * py[i];
        McDouble yf = -ck * px[i] + y[i] + sk * py[i];
        McDouble pyf = sn * px[i] + cs * py[i];

        // Update coords
        x[i] = xf;
        px[i] = pxf;
        y[i] = yf;
        py[i] = pyf;

        // Update the longitudinal position
        /*phifac = (pu_init * pu_init +
                  Bc * Bc * u_init * u_init +
                  dp_p * dp_p * gamma2i) / 2.0;
        phifac = (phifac * KNL - dp_p * gamma2i) * KNL;
        z[i] -= length * phifac;*/
    }
}

// Main tracking function
// Tracking function called by PyORBIT
void McMillan::trackBunch(Bunch* bunch, double length, double ks)
{
    int i, j; // Loop indices
    Nparticles = bunch->getSize(); // Number of particles
    McDouble slice_length = length/Nslices; // Length of each slice
    x.resize(Nparticles); y.resize(Nparticles); z.resize(Nparticles);
    px.resize(Nparticles); py.resize(Nparticles); pz.resize(Nparticles);
    flag.resize(Nparticles);

    McDouble k2 = 0.5*ks;

    // Step 1 and 2: Copy the bunch data into McDoubles and apply the entrance kick
    for (i = 0; i < Nparticles; i++)
    {
        // Copy the bunch data into internal McDouble vectors
        x[i] = bunch->x(i); y[i] = bunch->y(i); z[i] = bunch->z(i);
        px[i] = bunch->px(i); py[i] = bunch->py(i); pz[i] = bunch->pz(i);
        flag[i] = bunch->flag(i);
        // Apply the entrance rotation
        px[i] -= k2 * y[i];
        py[i] += k2 * x[i];
    }

    // Step 3: Loop over slices and perform symplectic integration steps
    for (i = 0; i < Nslices; i++)
    {
        // Loop through integrator
        for (j = 0; j < Ycoeffs->getNkicks(); j++)
        {
            PROPAGATEBUNCH(bunch, Ycoeffs->getTransportCoefficient(j)*slice_length, ks)
            kickBunch(Ycoeffs->getKickCoefficient(j), slice_length);
        }
        PROPAGATEBUNCH(bunch, Ycoeffs->getTransportCoefficient(j)*slice_length, ks)
    }

    // Step 4 and 5: Apply exit kick and copy the McDouble data back into the bunch
    for (i = 0; i < Nparticles; i++)
    {
        // Apply the exit rotation
        px[i] += k2 * y[i];
        py[i] -= k2 * x[i];
        // Copy the McDouble vectors back into the bunch
        bunch->x(i) = x[i]; bunch->y(i) = y[i]; bunch->z(i) = z[i];
        bunch->px(i) = px[i]; bunch->py(i) = py[i]; bunch->pz(i) = pz[i];
    }
}

/*****************************************************************************/
// Functions for the YoshidaCoefficients class

// Constructor - Supports all orders!
YoshidaCoefficients::YoshidaCoefficients(int Order)
{
    // Save parameters
    order = Order;
    Nkicks = std::pow(3, order/2 - 1);
    findCoefficients(); // construct coeffs
    
#ifdef VERBOSE 
    int i; double sum=0.0;
    std::cout << WARNING_PREFIX << "transport_coeffs[] = {";
    for(i=0;i<Nkicks+1;i++)
    {
        sum += transport_coeffs[i];
        std::cout << transport_coeffs[i] << ", ";
    }
    std::cout << "}, sum = " << sum << std::endl;

    sum = 0.0;
    std::cout << WARNING_PREFIX << "kick_coeffs[] = {";
    for(i=0;i<Nkicks;i++)
    {
        sum += kick_coeffs[i];
        std::cout << kick_coeffs[i] << ", ";
    }
    std::cout << "}, sum = " << sum << std::endl;
#endif
}


// Generalized construction of Yoshida Coefficients - valid for any order
void YoshidaCoefficients::findCoefficients()
{
    int i,j,k,steps = order / 2 - 1; // Number of iteration required to construct coeffs
    // Starting coefficients for order 2 - Most accurate choice.
    McDouble a = 0.5, b = 1.0, c = 0.5;
    // Final length of array containing all kick and transport coefficients.
    int full_L = 2 * std::pow(3, steps) + 1;
    McDouble K[full_L], Knew[full_L];

    // Initialize the seed values of order 2
    K[0] = a;
    K[1] = b;
    K[2] = c;

    // Interate through steps
    for (i = 0; i < steps; i++)
    {
        int L = 2 * std::pow(3, i) + 1; // Total number of coefficients for step i

        // Yoshida Magic!
        McDouble n = i + 1;
        McDouble z0 = -(std::pow(2, 1 / (2 * n + 1))/(2 - std::pow(2, 1 / (2 * n + 1))));
        McDouble z1 = (1 / (2 - std::pow(2, 1 / (2 * n + 1))));
        McDouble z;

        int c = 0; // Position in Knew
        for (j = 0; j < 3; j++) // Iterate over K three times - Magic!
        {
            if (j != 1) z = z1; // Alternate z among even and odd j
            else z = z0;

            for (k = 0; k < L; k++) // Iterate through the coefficients of the current order
                if (k != 0)
                {
                    Knew[c] = K[k] * z;
                    c += 1;
                }
                else if (c == 0)
                {
                    Knew[c] = K[k] * z;
                    c += 1;
                }
                else Knew[c-1] = Knew[c-1] + K[k] * z;
        }
        // Copy the new coefficients to the complete array
        for (j = 0; j < 3 * L - 2; j++) K[j] = Knew[j];
    }

    // Setting up the coefficient array for the gievn order
    int Num_kicks = std::pow(3, steps);
    kick_coeffs.resize(Num_kicks); // Resize vectors
    transport_coeffs.resize(Num_kicks + 1);

    for (int i = 0; i < Num_kicks; i++)
    {
        transport_coeffs[i] = K[2 * i];
        kick_coeffs[i] = K[2 * i + 1];
    }
    transport_coeffs[Num_kicks] = K[2 * Num_kicks];
}

/*****************************************************************************/
// Support functions
