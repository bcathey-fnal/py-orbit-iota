// Implements the tune buffer, see TuneBuffer.hh

#include "TuneBuffer.hh"

#include "orbit_mpi.hh"
#include "SyncPart.hh"
#include "ParticleAttributes.hh"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <string>

using namespace OrbitUtils;

static const double TWO_PI = 6.283185307179586;

// The class constants, defined here since they are bound by reference
const int TuneBuffer::NCOORD;
const int TuneBuffer::NSUMS;
const signed char TuneBuffer::UNKNOWN;

// Constructor
TuneBuffer::TuneBuffer(int window_in): CppPyWrapper(NULL)
{
    window = window_in;
    nslots = 0;
    latest = -1;
    turn = 0;
    maxstep = 0.0;
    double dnan = std::numeric_limits<double>::quiet_NaN();
    for(int p = 0; p < 2; p++)
    {
        cprevious[p] = cstart[p] = std::complex<double>(dnan, dnan);
        cphase[p] = 0.0;
    }
}

// Assign a slot to every particle of the bunch and size the arrays
void TuneBuffer::reset(Bunch* bunch, int window_in)
{
    int i;
    window = window_in;
    nslots = bunch->getSize();
    // The slot attribute, created if the bunch has none yet
    if(bunch->hasParticleAttributes("TuneSlot") == 0)
    {
        std::map<std::string, double> params;
        bunch->addParticleAttributes("TuneSlot", params);
    }
    ParticleAttributes* slot = bunch->getParticleAttributes("TuneSlot");
    ParticleAttributes* weight = NULL;
    if(bunch->hasParticleAttributes("macrosize") != 0)
        weight = bunch->getParticleAttributes("macrosize");
    isbeam.assign(nslots, 1);
    for(i = 0; i < nslots; i++)
    {
        slot->attValue(i, 0) = (double) i;
        if(weight != NULL && weight->attValue(i, 0) <= 0.0) isbeam[i] = 0;
    }
    float fnan = std::numeric_limits<float>::quiet_NaN();
    double dnan = std::numeric_limits<double>::quiet_NaN();
    coords.assign((size_t) NCOORD*nslots*window, fnan);
    windings.assign((size_t) 2*nslots*window, UNKNOWN);
    centroid.assign((size_t) NCOORD*window, dnan);
    cwindings.assign((size_t) 2*window, UNKNOWN);
    turns.assign(window, -1L);
    previous.assign((size_t) 2*nslots, std::complex<float>(fnan, fnan));
    start.assign((size_t) 2*nslots, std::complex<float>(fnan, fnan));
    phase.assign((size_t) 2*nslots, 0.0f);
    for(int p = 0; p < 2; p++)
    {
        cprevious[p] = cstart[p] = std::complex<double>(dnan, dnan);
        cphase[p] = 0.0;
    }
    recorded.clear();
    latest = -1;
    turn = 0;
    maxstep = 0.0;
}

// Moments of the beam particles where the bunch is, reduced over MPI
void TuneBuffer::moments(Bunch* bunch, double* sums)
{
    int i, a, b, n = bunch->getSize();
    SyncPart* sp = bunch->getSyncPart();
    double dpfac = 1.0/(sp->getMomentum()*sp->getBeta()); // dE to delta
    ParticleAttributes* weight = NULL;
    if(bunch->hasParticleAttributes("macrosize") != 0)
        weight = bunch->getParticleAttributes("macrosize");
    double local[NSUMS], u[NCOORD];
    for(a = 0; a < NSUMS; a++) local[a] = 0.0;
    for(i = 0; i < n; i++)
    {
        if(weight != NULL && weight->attValue(i, 0) <= 0.0) continue;
        u[0] = bunch->x(i); u[1] = bunch->xp(i);
        u[2] = bunch->y(i); u[3] = bunch->yp(i);
        u[4] = bunch->z(i); u[5] = bunch->dE(i)*dpfac;
        local[0] += 1.0;
        for(a = 0; a < NCOORD; a++)
        {
            local[1 + a] += u[a];
            for(b = 0; b < NCOORD; b++)
                local[1 + NCOORD + NCOORD*a + b] += u[a]*u[b];
        }
    }
    ORBIT_MPI_Allreduce(local, sums, NSUMS, MPI_DOUBLE, MPI_SUM,
                        bunch->getMPI_Comm_Local()->comm);
}

// The closed-form normalization from the moments
bool TuneBuffer::normalization(const double* sums, double* mean,
                               std::complex<double>* cvec)
{
    int a, b, i, j, c, l;
    double dnan = std::numeric_limits<double>::quiet_NaN();
    for(c = 0; c < 2*NCOORD; c++) cvec[c] = std::complex<double>(dnan, dnan);
    double count = sums[0];
    if(count < 2.0) return false;
    double sigma[NCOORD][NCOORD];
    for(a = 0; a < NCOORD; a++) mean[a] = sums[1 + a]/count;
    for(a = 0; a < NCOORD; a++)
        for(b = 0; b < NCOORD; b++)
            sigma[a][b] = sums[1 + NCOORD + NCOORD*a + b]/count -
                          mean[a]*mean[b];
    for(a = 0; a < NCOORD; a++)
        for(b = 0; b < a; b++)
        {
            double s = 0.5*(sigma[a][b] + sigma[b][a]);
            sigma[a][b] = sigma[b][a] = s;
        }
    // The regression of (x, x', y, y') on (z, delta)
    double S00 = sigma[4][4], S01 = sigma[4][5], S11 = sigma[5][5];
    double det = S00*S11 - S01*S01;
    double D[4][2];
    for(i = 0; i < 4; i++) D[i][0] = D[i][1] = 0.0;
    if(det > 0.0 && S00 > 0.0 && S11 > 0.0)
    {
        double inv00 = S11/det, inv01 = -S01/det, inv11 = S00/det;
        for(i = 0; i < 4; i++)
        {
            D[i][0] = sigma[i][4]*inv00 + sigma[i][5]*inv01;
            D[i][1] = sigma[i][4]*inv01 + sigma[i][5]*inv11;
        }
    }
    // The betatron block, the correlated part removed
    double block[4][4];
    for(i = 0; i < 4; i++)
        for(j = 0; j < 4; j++)
            block[i][j] = sigma[i][j] -
                (D[i][0]*(S00*D[j][0] + S01*D[j][1]) +
                 D[i][1]*(S01*D[j][0] + S11*D[j][1]));
    bool ok = true;
    for(j = 0; j < 2; j++)
    {
        i = 2*j;
        double s11 = block[i][i], s12 = block[i][i + 1],
               s22 = block[i + 1][i + 1];
        double eps2 = s11*s22 - s12*s12;
        if(eps2 <= 0.0 || s11 <= 0.0)
        {
            ok = false;
            continue;
        }
        double eps = std::sqrt(eps2);
        double beta = s11/eps, alpha = -s12/eps, sb = std::sqrt(beta);
        double T0[NCOORD], T1[NCOORD];
        for(c = 0; c < NCOORD; c++) T0[c] = T1[c] = 0.0;
        T0[i] = 1.0/sb;
        T1[i] = alpha/sb;
        T1[i + 1] = sb;
        for(l = 0; l < 2; l++)
        {
            T0[4 + l] = -(T0[i]*D[i][l] + T0[i + 1]*D[i + 1][l]);
            T1[4 + l] = -(T1[i]*D[i][l] + T1[i + 1]*D[i + 1][l]);
        }
        for(c = 0; c < NCOORD; c++)
            cvec[j*NCOORD + c] = std::complex<double>(T0[c], -T1[c]);
    }
    return ok;
}

// The tune node's action
void TuneBuffer::trackBunch(Bunch* bunch, bool isfirst)
{
    if(nslots == 0) return; // reset has not been called
    int i, c, p, n = bunch->getSize();
    if(bunch->hasParticleAttributes("TuneSlot") == 0)
        ORBIT_MPI_Finalize("monitor.tunebuffer: the bunch has no TuneSlot"
                           " attribute; call reset(bunch, window) first.");
    ParticleAttributes* slotatt = bunch->getParticleAttributes("TuneSlot");
    SyncPart* sp = bunch->getSyncPart();
    double dpfac = 1.0/(sp->getMomentum()*sp->getBeta());
    double sums[NSUMS], mean[NCOORD], u[NCOORD];
    std::complex<double> cvec[2*NCOORD];
    moments(bunch, sums);
    bool ok = normalization(sums, mean, cvec);
    float fnan = std::numeric_limits<float>::quiet_NaN();
    double dnan = std::numeric_limits<double>::quiet_NaN();
    int col = -1, prevcol = -1;
    if(isfirst)
    {
        col = (int) (turn % window);
        // The column starts empty: the particles present fill it
        for(c = 0; c < NCOORD; c++)
            for(i = 0; i < nslots; i++)
                coords[((size_t) c*nslots + i)*window + col] = fnan;
        for(p = 0; p < 2; p++)
            for(i = 0; i < nslots; i++)
                windings[((size_t) p*nslots + i)*window + col] = UNKNOWN;
        if(turn > 0) prevcol = (int) ((turn - 1) % window);
        recorded.assign(n, -1);
    }
    for(i = 0; i < n; i++)
    {
        int s = (int) slotatt->attValue(i, 0);
        if(s < 0 || s >= nslots) continue; // Not a particle of the reset
        u[0] = bunch->x(i); u[1] = bunch->xp(i);
        u[2] = bunch->y(i); u[3] = bunch->yp(i);
        u[4] = bunch->z(i); u[5] = bunch->dE(i)*dpfac;
        if(isfirst)
        {
            for(c = 0; c < NCOORD; c++)
                coords[((size_t) c*nslots + s)*window + col] = (float) u[c];
            recorded[i] = s;
        }
        if(!ok) continue; // No normalization here: no step
        for(p = 0; p < 2; p++)
        {
            std::complex<double> hd(0.0, 0.0);
            for(c = 0; c < NCOORD; c++)
                hd += cvec[p*NCOORD + c]*(u[c] - mean[c]);
            std::complex<float> h((float) hd.real(), (float) hd.imag());
            size_t k = (size_t) p*nslots + s;
            std::complex<float> prev = previous[k];
            if(!std::isnan(prev.real())) // Seen before: add the step
            {
                double d = std::arg(std::complex<double>(h)*
                                    std::conj(std::complex<double>(prev)));
                phase[k] += (float) d;
                if(std::fabs(d) > maxstep) maxstep = std::fabs(d);
            }
            previous[k] = h;
            if(isfirst)
            {
                // Close the turn: the phase accumulated over it less the
                // direct start-to-start step is a whole number of 2 pi
                if(prevcol >= 0 && !std::isnan(start[k].real()))
                {
                    double direct = std::arg(std::complex<double>(h)*
                                    std::conj(std::complex<double>(start[k])));
                    double kt = std::floor((phase[k] - direct)/TWO_PI + 0.5);
                    if(kt < -127.0) kt = -127.0;
                    if(kt > 127.0) kt = 127.0;
                    windings[k*window + prevcol] = (signed char) kt;
                }
                phase[k] = 0.0f;
                start[k] = h;
            }
        }
    }
    // The centroid's chain, about the design orbit
    if(ok)
    {
        for(p = 0; p < 2; p++)
        {
            std::complex<double> hc(0.0, 0.0);
            for(c = 0; c < NCOORD; c++) hc += cvec[p*NCOORD + c]*mean[c];
            if(!std::isnan(cprevious[p].real()))
                cphase[p] += std::arg(hc*std::conj(cprevious[p]));
            cprevious[p] = hc;
            if(isfirst)
            {
                if(prevcol >= 0 && !std::isnan(cstart[p].real()))
                {
                    double direct = std::arg(hc*std::conj(cstart[p]));
                    double kt = std::floor((cphase[p] - direct)/TWO_PI + 0.5);
                    if(kt < -127.0) kt = -127.0;
                    if(kt > 127.0) kt = 127.0;
                    cwindings[(size_t) p*window + prevcol] = (signed char) kt;
                }
                cphase[p] = 0.0;
                cstart[p] = hc;
            }
        }
    }
    if(isfirst)
    {
        for(c = 0; c < NCOORD; c++)
            centroid[(size_t) c*window + col] = ok ? mean[c] : dnan;
        for(p = 0; p < 2; p++) cwindings[(size_t) p*window + col] = UNKNOWN;
        turns[col] = turn;
        latest = col;
        turn++;
    }
}
