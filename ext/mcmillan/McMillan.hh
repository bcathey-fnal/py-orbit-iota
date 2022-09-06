// Declares the C++ McMillan thick lens class and associated functions

#ifndef MCMILLAN_H
#define MCMILLAN_H

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
#include <cmath>

// Useful macros
#define WARNING_PREFIX "Warning from line " << __LINE__ << " in file " << __FILE__ << ": "
#define PROPAGATEBUNCH(B,L,KS) if(abs(KS) < KSMIN) driftBunch((B), (L)); else solenoidBunch((B), (L), (KS));

// Floating point type for McMillan lens calculation
typedef double McDouble; // Too funny not to do it! ><

// Hard coded parameters - Add any algorithm or other physical parameters as macros
// eg. #define kB 1.38064852e-23 // Boltzmann constant in J/K
#define KSMIN 1e-4 // Minimum value of ks below which the solenoid is a drift

/*****************************************************************************/
// Class for Symplectic Integration coefficients
class YoshidaCoefficients
{
    private:
        int order, Nkicks;
        std::vector<McDouble> kick_coeffs, transport_coeffs;
        void findCoefficients(); // Fill in coefficient vectors

    public:
        // Constructor for coefficients
        YoshidaCoefficients(int Order);
        ~YoshidaCoefficients(){}; // Destructor
        inline int getOrder(){return order;}; // Get the order
        inline int getNkicks(){return Nkicks;}; // Get the number of kicks
        // Get kick coefficient
        inline McDouble getKickCoefficient(int index){return kick_coeffs[index];};
        // Get drift coefficient
        inline McDouble getTransportCoefficient(int index)
        {
            return transport_coeffs[index];
        };
};

/*****************************************************************************/
// McMillan lens class
class McMillan: public OrbitUtils::CppPyWrapper
{
    private:
        // All parameters for the McMillan lens
        double strength_over_length, width;
        int Nslices, Nparticles;
        // Store the symplectic coefficient object
        YoshidaCoefficients *Ycoeffs;

        // Local storage of bunch data
        std::vector<McDouble> x,y,z,px,py,pz;
        std::vector<int> flag;

        // McMillan kicker
        void kickBunch(McDouble coeff, McDouble slice_length);
        // Copy of the drift teapot code
        void driftBunch(Bunch* bunch, McDouble length);
        // Drift in the presence of solenoid field
        void solenoidBunch(Bunch* bunch, McDouble length, McDouble ks);


    public:
        // Constructor. K: strength_over_length, A: width, O:order, N: slices
        McMillan(double K, double A, int O, int N);
        // Destructor - Don't worry about this yet.
        ~McMillan();
        // Basic functions to access class members
        inline double getStrengthOverLength(){return strength_over_length;};
        inline double getWidth(){return width;};
        inline int getNslices(){return Nslices;};
        inline int getOrder(){return Ycoeffs->getOrder();};
        void setNslices(int N){Nslices=N;};

        // Track the bunch passing through the lens - Called from the
        // McMillan_AccNode class
        void trackBunch(Bunch* bunch, double length, double ks);
};

/*****************************************************************************/
// Additional support functions

#endif // End of header

