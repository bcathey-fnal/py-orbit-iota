//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   TuneSlot.hh
//
// DESCRIPTION
//    A subclass of the particle attributes class: the slot of a
//    macro-particle in the ring buffer of the tune diagnostics
//    (ext/monitor/TuneBuffer). Assigned once, it follows the particle through
//    compressions of the bunch and into the lost bunch, which is what lets
//    turn-by-turn records be attributed to particles after losses.
//
///////////////////////////////////////////////////////////////////////////
#ifndef TUNE_SLOT_H
#define TUNE_SLOT_H

#include <string>

#include "ParticleAttributes.hh"

class TuneSlot : public ParticleAttributes
{
public:
  TuneSlot(Bunch* bunch);
  ~TuneSlot();
  /** Returns the slot of the particle with index. */
  int getSlot(int particle_index);
  /** Sets the slot of the particle with index. */
  void setSlot(int particle_index, int slot);
};

#endif
