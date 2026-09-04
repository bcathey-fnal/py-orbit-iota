//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    OrbitConst.cc
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    06/28/2005
//
// DESCRIPTION
//    Keeps all constants for ORBIT
//
///////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//
// include files
//
/////////////////////////////////////////////////////////////////////////////


#include "OrbitConst.hh"

const double OrbitConst::PI = 3.14159265358979323846264;
const double OrbitConst::c = 2.99792458e+8;
const double OrbitConst::elementary_charge_CGS = 4.803204712570263e-10;
const double OrbitConst::coeff_Tesla_to_inner = (1.0e+4/4.803204712570263e-10)*0.01*(1.0e+6);
const double OrbitConst::coeff_VoltsPerM_to_inner = 1.0/((1.602176634e-19)*(8.987551786170797e+9));
const double OrbitConst::coeff_Phi_to_Volts = (1.602176634e-19)*(8.987551786170797e+9);
const double OrbitConst::mass_electron = 0.51099895069e-3;
const double OrbitConst::classicalRadius_electron = 2.8179403205e-15;
const double OrbitConst::charge_electron = -1.0;
const double OrbitConst::mass_proton = 0.93827208943;
const double OrbitConst::classicalRadius_proton = 1.5346982641e-18;
const double OrbitConst::charge_proton = +1.0;
const double OrbitConst::permeability = 1.25663706127e-6;
const double OrbitConst::elementary_charge_MKS = 1.602176634e-19;
const double OrbitConst::tiny = 2.22044604925031e-16;


OrbitConst::OrbitConst()
{
}

OrbitConst::~OrbitConst()
{
}

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////
