// Implements the beam distribution monitor
#ifndef BEAMDIST_H
#define BEAMDIST_H

// ORBIT includes
#include "Bunch.hh"
#include "OrbitConst.hh"
// pyORBIT utils
#include "CppPyWrapper.hh"

// Standard includes
#include <iostream>
#include <vector>

// Define some contants to hold trigger type
typedef enum {TRIG_SINGLE, TRIG_RUN} TRIG_TYPE;

/*****************************************************************************/
// Beam distribution class
class BeamDist: public OrbitUtils::CppPyWrapper
{
    private:
        TRIG_TYPE trigger_mode; // Which type of trigger to use
        int nparticles, trigger_line; // Number of particles and trigger line
        std::vector<double> x, y, z, xp, yp, dE; // Particle coordinates

    public:
        BeamDist(int npart, TRIG_TYPE trig_mode); // Constructor
        ~BeamDist(); // Destructor
        void trigger(); // Trigger an acquisition
        // Track the bunch passing through the monitor
        void trackBunch(Bunch* bunch);
        // Get references to the particle data
        std::vector<double>& getx() {return x;};
        std::vector<double>& gety() {return y;};
        std::vector<double>& getz() {return z;};
        std::vector<double>& getxp() {return xp;};
        std::vector<double>& getyp() {return yp;};
        std::vector<double>& getdE() {return dE;};
};

#endif // End of header
