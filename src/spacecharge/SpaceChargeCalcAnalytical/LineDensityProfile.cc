/////////////////////////////////////////////////////////////////////////////
//
// FILE NAME
//   LineDensityProfile.cc
//
//   17/02/2015
//
// DESCRIPTION
//   Provides Linedensity Profiles needed for the analytical SC calculator
//
/////////////////////////////////////////////////////////////////////////////

#include "LineDensityProfile.hh"
#include "BufferStore.hh"
#include "OrbitConst.hh"
#include "CppPyWrapper.hh"

#include <iostream>
#include <cmath>
#include <cfloat>
#include <complex>

using namespace OrbitUtils;


// Class of LineDensity Profiles
LineDensityProfile::LineDensityProfile() : CppPyWrapper(NULL) {};

double LineDensityProfile::getLocalLineDensityFactor(double z) 
{
	return 1.;
};

// const char* LineDensityProfile::getType() 
// {
// 	return type
// };

// Derived Class for Gaussian LineDensity
GaussianLineDensityProfile::GaussianLineDensityProfile(double blength_rms): LineDensityProfile() 
{
	blength_rms_ = blength_rms;
// 	int rank = 0;
// 	ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
// 	if(rank == 0){std::cerr << "  SpaceChargeCalcAnalytical::GaussianLineDensityProfile initialized (rms bunch length = " << blength_rms_ << "m)" << std::endl;};
};

double GaussianLineDensityProfile::getLocalLineDensityFactor(double z) 
{
	return 1. / ( blength_rms_ * sqrt(2*OrbitConst::PI) ) * exp(-pow(z, 2) / 2 / pow(blength_rms_, 2) );
};

void GaussianLineDensityProfile::setBunchLength(double blength_rms)
{
	blength_rms_ = blength_rms;
};


// Derived Class for Constant LineDensity
ConstantLineDensityProfile::ConstantLineDensityProfile(double blength_full): LineDensityProfile() 
{
	blength_full_ = blength_full;
// 	int rank = 0;
// 	ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
// 	if(rank == 0){std::cerr << "  SpaceChargeCalcAnalytical::ConstantLineDensityProfile initialized (full bunch length = " << blength_full_ << "m)" << std::endl;};
};

double ConstantLineDensityProfile::getLocalLineDensityFactor(double z) 
{
	return 1. / blength_full_;
};

void ConstantLineDensityProfile::setBunchLength(double blength_full)
{
	blength_full_ = blength_full;
};


// Derived Class for Interpolated LineDensity
InterpolatedLineDensityProfile::InterpolatedLineDensityProfile(double z_min, double z_max, double lambda[], int array_length): LineDensityProfile() 
{
	this->init(z_min, z_max, lambda, array_length);
};

// Destructor
InterpolatedLineDensityProfile::~InterpolatedLineDensityProfile()
{
	delete [] z_;
	delete [] lambda_;
}


void InterpolatedLineDensityProfile::init(double z_min, double z_max, double lambda[], int array_length)
{
	if ( !(z_max > z_min) ){ORBIT_MPI_Finalize("PySpaceChargeCalcAnalytical.setInterpolatedLineDensityProfile(z_min, z_max, lambda[]) - z_max > z_min required.");}
	if (array_length < 2) {ORBIT_MPI_Finalize("Error! SpaceChargeCalcAnalytical::InterpolateLineDensityProfile expects arrays of length > 1 ");}
		
	z_min_ = z_min;
	z_max_ = z_max;
	dz_ = (z_max - z_min) / (array_length - 1);
	array_length_ = array_length;
	
	z_ = new double[array_length_];
	lambda_ = new double[array_length_];
	
	for (int i=0; i<array_length_; i++) {
		z_[i] = z_min_ + i * dz_;
		lambda_[i] = lambda[i];
		//std::cerr << "i = " << i << ", z_[i] = " << z_[i] << ", lambda_[i] = " << lambda_[i] << std::endl;
	};	
}


void InterpolatedLineDensityProfile::setLineDensityProfile(double z_min, double z_max, double lambda[], int array_length)
{
	delete [] z_;
	delete [] lambda_;
	this->init(z_min, z_max, lambda, array_length);
}


double InterpolatedLineDensityProfile::getLocalLineDensityFactor(double z) 
{
	if (z < z_min_ | z > z_max_) { return 0.;}
	else {
		double fractional_index = (z - z_min_) / dz_;
		int index = int(fractional_index);
		fractional_index -= index;
		//std::cerr << " index = " << index << ", fractional_index = " << fractional_index << std::endl;
		return (1 - fractional_index) * lambda_[index] + fractional_index * lambda_[index+1];
	}
};

