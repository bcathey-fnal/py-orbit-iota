// The ring buffer and the winding counter of the tune diagnostics of the
// IOTA space-charge benchmark: the C++ side of what Simulations/tunespectra.py
// of that benchmark asks a code to record, done where PyORBIT's bunch lives.
//
// One instance serves every tune node of the lattice. The node at the start
// of the ring calls trackBunch with isfirst = true, which closes the count of
// whole betatron oscillations of the turn just tracked and records the bunch
// into the buffer; every other node adds each particle's phase advance since
// the previous position. At each position the coordinates are normalized
// with the beam's own moments there, in closed form: the longitudinal pair
// (z, delta) regressed out of the transverse coordinates, then beta and alpha
// of each plane from its 2x2 betatron block. The moments are reduced over
// MPI, so every rank must call trackBunch at every node.
//
// The counter need not run on every turn. With setClose(false) the node at
// the start of the ring restarts the count instead of closing it: the turn
// just tracked keeps unknown windings and the step across it is not added,
// so a caller that switches the other nodes off for some turns sets close to
// whether the turn just tracked was counted.
//
// Particles are addressed by the TuneSlot particle attribute, which reset
// assigns from the bunch index at that moment and which follows a particle
// through compressions of the bunch. The arrays are laid out as
// Simulations/tunespectra.py's RingBuffer holds them, and python views them
// without copying.

#ifndef TUNEBUFFER_H
#define TUNEBUFFER_H

// ORBIT includes
#include "Bunch.hh"
// pyORBIT utils
#include "CppPyWrapper.hh"

// Standard includes
#include <complex>
#include <vector>

class TuneBuffer: public OrbitUtils::CppPyWrapper
{
    public:
        static const int NCOORD = 6; // x, x', y, y', z, delta
        static const int NSUMS = 1 + NCOORD + NCOORD*NCOORD; // Moments
        static const signed char UNKNOWN = -128; // A winding not counted

        TuneBuffer(int window); // Constructor
        ~TuneBuffer(){}; // Destructor
        // Assign a slot to every particle of the bunch and size the arrays
        void reset(Bunch* bunch, int window);
        // The tune node's action: count, and record at the start of the ring
        void trackBunch(Bunch* bunch, bool isfirst);
        // Whether the node at the start of the ring closes the turn just
        // tracked (true, the default) or restarts the count (false)
        void setClose(bool close_in) {close = close_in;};
        bool getClose() {return close;};

        // Sizes and counters
        int getNSlots() {return nslots;};
        int getWindow() {return window;};
        long getTurn() {return turn;}; // Turns recorded so far
        int getLatest() {return latest;}; // Column of the latest turn
        double getLargestStep() {return maxstep;}; // Largest |phase step|

        // The arrays, which python views as numpy arrays
        // (6, nslots, window): the coordinates, NaN where not recorded
        std::vector<float>& getCoords() {return coords;};
        // (2, nslots, window): whole windings over the turn of a column
        std::vector<signed char>& getWindings() {return windings;};
        // (6, window): the centroid of the beam on each turn
        std::vector<double>& getCentroid() {return centroid;};
        // (2, window): whole windings of the centroid
        std::vector<signed char>& getCentroidWindings() {return cwindings;};
        // (window,): the turn each column holds, -1 when empty
        std::vector<long>& getTurns() {return turns;};
        // (nslots,): 1 for beam particles, 0 for test particles
        std::vector<unsigned char>& getIsBeam() {return isbeam;};
        // The slots of the particles in bunch order at the last record
        std::vector<int>& getRecordedSlots() {return recorded;};

    private:
        int window, nslots, latest;
        long turn;
        bool close; // The start of the ring closes the turn, or restarts
        std::vector<float> coords;
        std::vector<signed char> windings;
        std::vector<double> centroid;
        std::vector<signed char> cwindings;
        std::vector<long> turns;
        std::vector<unsigned char> isbeam;
        std::vector<int> recorded;
        // The counter's state, index plane*nslots + slot: h at the previous
        // position and at the start of the ring, and the phase accumulated
        // since the start
        std::vector<std::complex<float> > previous, start;
        std::vector<float> phase;
        // The same for the centroid, about the design orbit
        std::complex<double> cprevious[2], cstart[2];
        double cphase[2];
        double maxstep;

        // Moments of the beam particles where the bunch is, reduced over MPI
        void moments(Bunch* bunch, double* sums);
        // The closed-form normalization from the moments: the mean and the
        // two complex vectors h_p = cvec[p] . (u - mean). False when the
        // beam has no emittance in a plane.
        bool normalization(const double* sums, double* mean,
                           std::complex<double>* cvec);
};

#endif // End of header
