#ifndef WRAP_EXT_H
#define WRAP_EXT_H

#include "wrap_monitor.hh"
#include "wrap_mcmillan.hh"
#include "wrap_scattering.hh"

// Register all extension modules
void initext()
{
    wrap_monitor::initmonitor();
    wrap_mcmillan::initmcmillan();
    wrap_scattering::initscattering();
}

#endif // WRAP_EXT_H

