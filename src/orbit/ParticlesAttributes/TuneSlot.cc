//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   TuneSlot.cc
//
// DESCRIPTION
//    A subclass of a ParticleAttributes class to keep the slot of each
//    particle in the ring buffer of the tune diagnostics.
//
///////////////////////////////////////////////////////////////////////////
#include "Bunch.hh"
#include "TuneSlot.hh"

TuneSlot::TuneSlot(Bunch* bunch):
  ParticleAttributes(bunch,1)
{
  cl_name_ = "TuneSlot";
  attrDescr = "Tune_Slot";
}

TuneSlot::~TuneSlot()
{
}

/** Returns the slot of the particle with index. */
int TuneSlot::getSlot(int particle_index)
{
	return int(attValue(particle_index,0));
}

/** Sets the slot of the particle with index. */
void TuneSlot::setSlot(int particle_index, int slot)
{
	attValue(particle_index,0) = 1.0*slot;
}
