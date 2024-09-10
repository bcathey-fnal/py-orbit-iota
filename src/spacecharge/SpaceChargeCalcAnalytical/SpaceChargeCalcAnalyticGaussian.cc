/////////////////////////////////////////////////////////////////////////////
//
// FILE NAME
//   SpaceChargeCalcAnalyticGaussian.cc
//
// AUTHOR
//   H. Bartosik
//
//   02/20/2015
//
// DESCRIPTION
//   Calculate the space charge effect of the bunch using a frozen space 
//   charge potential for a transverse bi-Gaussian distribution. For more
//   details see M. Bassetti and G.A. Erskine, "Closed expression for the
//   electric field of a two-dimensional Gaussian charge density", 
//   CERN-ISR-TH/80-06.
//
/////////////////////////////////////////////////////////////////////////////

#include "SpaceChargeCalcAnalyticGaussian.hh"
#include "BufferStore.hh"
#include "OrbitConst.hh"
#include "LineDensityProfile.hh"
#include "CppPyWrapper.hh"
#include "ErrorFunctions.hh"

#include <iostream>
#include <cmath>
#include <cfloat>
#include <complex>
#include <typeinfo>

using namespace OrbitUtils;

SpaceChargeCalcAnalyticGaussian::SpaceChargeCalcAnalyticGaussian(double intensity, double epsn_x, double epsn_y, double dpp_rms, LineDensityProfile* LineDensityProfile): CppPyWrapper(NULL)
{
	LineDensityProfile_ = NULL;
	this->setBunchParameters(intensity, epsn_x, epsn_y, dpp_rms);
	this->setLineDensityProfile(LineDensityProfile);
}


SpaceChargeCalcAnalyticGaussian::~SpaceChargeCalcAnalyticGaussian(){
	Py_DECREF(LineDensityProfile_->getPyWrapper());
};


void SpaceChargeCalcAnalyticGaussian::setLineDensityProfile(LineDensityProfile* LineDensityProfile)
{
	if (LineDensityProfile_ != NULL){
		Py_DECREF(LineDensityProfile_->getPyWrapper());
	};
	LineDensityProfile_ = LineDensityProfile;
	Py_INCREF(LineDensityProfile_->getPyWrapper());
};


void SpaceChargeCalcAnalyticGaussian::setBunchParameters(double intensity, double epsn_x, double epsn_y, double dpp_rms){
	intensity_ = intensity;
	epsn_x_ = epsn_x;
	epsn_y_ = epsn_y;	
	dpp_rms_ = dpp_rms;
}


void SpaceChargeCalcAnalyticGaussian::trackBunch(Bunch* bunch, double length){
	SyncPart* syncPart = bunch->getSyncPart();
	double beta = syncPart->getBeta();
	double gamma = syncPart->getGamma(); 
		
	double sigma_x = sqrt(epsn_x_ * beta_x_ / beta / gamma + dpp_rms_*dpp_rms_ * eta_x_*eta_x_);
	double sigma_y = sqrt(epsn_y_ * beta_y_ / beta / gamma + dpp_rms_*dpp_rms_ * eta_y_*eta_y_);
	
	double SC_factor = length * intensity_ * 4 * OrbitConst::PI * bunch->getClassicalRadius() * pow(bunch->getCharge(),2) / (pow(beta, 2) * pow(gamma, 3));

	double x, y, z, fx, fy, L_factor;
	for (int i = 0, n = bunch->getSize(); i < n; i++){
		x = bunch->x(i) - co_x_;
		y = bunch->y(i) - co_y_;
		z = bunch->z(i);
		this->BassettiErskine(x, y, sigma_x, sigma_y, fx, fy);
		L_factor = this->getLineDensityFactor(z);
		bunch->xp(i) += fx * SC_factor * L_factor;
		bunch->yp(i) += fy * SC_factor * L_factor;	
	}		
}


void SpaceChargeCalcAnalyticGaussian::setLatticeParameters(double beta_x, double beta_y, double eta_x, double eta_y, double co_x, double co_y){
	beta_x_ = beta_x;
	beta_y_ = beta_y;
	eta_x_  = eta_x;
	eta_y_  = eta_y;
	co_x_   = co_x;
	co_y_   = co_y;		
}


double SpaceChargeCalcAnalyticGaussian::getLineDensityFactor(double z) 
{
	return LineDensityProfile_->getLocalLineDensityFactor(z);
};


void SpaceChargeCalcAnalyticGaussian::BassettiErskine(double x, double y, double sigma_x, double sigma_y, double& E_x, double& E_y){
	double eps0=1.;
	int sign_x, sign_y;
	
	if (x < 0){ sign_x = -1;  x = abs(x); }
	else{ sign_x = 1; }
	
	if (y < 0){ sign_y = -1;  y = abs(y); }
	else{ sign_y = 1; }

	double S, factor_BE;
	if (sigma_x > sigma_y){
		S = sqrt(2 * (sigma_x * sigma_x - sigma_y * sigma_y));
		factor_BE = 1 / (2 * eps0 * sqrt(OrbitConst::PI) * S);
		std::complex<double> arg1(sigma_y / sigma_x * x, sigma_x / sigma_y * y);
		std::complex<double> arg2(x, y);

		std::complex<double> E = factor_BE * (cerrf(arg2/S) - exp( -x * x / (2 * sigma_x * sigma_x) - y * y / (2 * sigma_y * sigma_y)) * cerrf(arg1/S) );

		E_x = abs(imag(E)) * sign_x;
		E_y = abs(real(E)) * sign_y;
	}
	else if (sigma_x < sigma_y){
		S = sqrt(2 * (sigma_y * sigma_y - sigma_x * sigma_x));
		factor_BE = 1 / (2 * eps0 * sqrt(OrbitConst::PI) * S);
		std::complex<double> arg1(sigma_x / sigma_y * y, sigma_y / sigma_x * x);
		std::complex<double> arg2(y, x);
		
		std::complex<double> E = factor_BE * (cerrf(arg2/S) - exp( -y * y / (2 * sigma_y * sigma_y) - x * x / (2 * sigma_x * sigma_x)) * cerrf(arg1/S) );
		
		E_y = abs(imag(E)) * sign_y;
		E_x = abs(real(E)) * sign_x;
	}
	else {
		if (x == 0 & y ==0) {
			E_x = 0.;
			E_y = 0.;
		} else {
			double r2 = pow(x, 2) + pow(y, 2);
			factor_BE = 1 / (2 * OrbitConst::PI * eps0 * r2) * (1 - exp(- r2 / (pow(sigma_x,2) + pow(sigma_y, 2)) ) );
			E_x = factor_BE * x * sign_x;
			E_y = factor_BE * y * sign_y;			
		}
	}
}