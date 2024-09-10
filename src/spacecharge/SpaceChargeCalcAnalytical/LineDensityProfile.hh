//   Class for Linedensity Profiles needed for the analytical SC calculator

#ifndef LINEDENSITY_PROFILES_H
#define LINEDENSITY_PROFILES_H

//MPI Function Wrappers
#include "orbit_mpi.hh"
#include "wrap_mpi_comm.hh"

#include <cstdlib>
#include <cmath>

//ORBIT bunch
#include "Bunch.hh"
//pyORBIT utils
#include "CppPyWrapper.hh"

using namespace std;



class LineDensityProfile: public OrbitUtils::CppPyWrapper
{ 
	public:
		//* Constructor *
	 	LineDensityProfile();
 	
		virtual double getLocalLineDensityFactor(double z);
};
 
class GaussianLineDensityProfile : public LineDensityProfile 
{
	public:
		//* Constructor *
		GaussianLineDensityProfile(double blength_rms);
 
		double getLocalLineDensityFactor(double z);
		
		void setBunchLength(double blength_rms);
 		
	protected: 
		double blength_rms_; 		
};


class ConstantLineDensityProfile : public LineDensityProfile 
{
	public:
		//* Constructor *
 		ConstantLineDensityProfile(double blength_full);
 		
		double getLocalLineDensityFactor(double z);
		
		void setBunchLength(double blength_full);
 		
 	protected: 
		double blength_full_; 		
};
	

class InterpolatedLineDensityProfile : public LineDensityProfile 
{
	public:
		//* Constructor *
		InterpolatedLineDensityProfile(double z_min, double z_max, double lambda[], int array_length);
 		
 		//* Destructor *
 		~InterpolatedLineDensityProfile();
 		
 		void init(double z_min, double z_max, double lambda[], int array_length);
 		
 		void setLineDensityProfile(double z_min, double z_max, double lambda[], int array_length);
 		
		double getLocalLineDensityFactor(double z);
 		
 	protected: 
		double z_min_;
		double z_max_;
 		double dz_;
 		int array_length_;
		double* z_;
		double* lambda_;
};
#endif
