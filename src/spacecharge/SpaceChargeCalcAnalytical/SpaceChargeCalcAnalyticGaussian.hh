//   Calculate the space charge effect of the bunch using a frozen space charge potential  

#ifndef SC_SPACECHARGE_ANALYTIC_GAUSSIAN_H
#define SC_SPACECHARGE_ANALYTIC_GAUSSIAN_H

//MPI Function Wrappers
#include "orbit_mpi.hh"
#include "wrap_mpi_comm.hh"

#include <cstdlib>
#include <cmath>

//ORBIT bunch
#include "Bunch.hh"
//pyORBIT utils
#include "CppPyWrapper.hh"

#include "LineDensityProfile.hh"

using namespace std;


class SpaceChargeCalcAnalyticGaussian: public OrbitUtils::CppPyWrapper {
	
	public:
	
		// Constructor
		SpaceChargeCalcAnalyticGaussian(double intensity, double epsn_x, double epsn_y, double dpp_rms, LineDensityProfile* LineDensityProfile);
	
		// Destructor
		~SpaceChargeCalcAnalyticGaussian();
		
		// for setting the bunch parameters
		void setBunchParameters(double intensity, double epsn_x, double epsn_y, double dpp_rms);
	
		// to update the line density profile object
		void setLineDensityProfile(LineDensityProfile* LineDensityProfile);
			
		/** Calculates space charge and applies the transverse SC kicks 
		to the macro-particles in the bunch. */
		void trackBunch(Bunch* bunch, double length);

		void setLatticeParameters(double beta_x, double beta_y, double eta_x, double eta_y, double co_x, double co_y);
		
		double getLineDensityFactor(double z);
				
		void BassettiErskine(double x, double y, double sigma_x, double sigma_y, double& E_x, double& E_y);
		
	private:
		 		

	protected:
	
		double intensity_;
		double epsn_x_;
		double epsn_y_;
		double dpp_rms_;
		double beta_x_;
		double beta_y_;
		double eta_x_;
		double eta_y_;
		double co_x_;
		double co_y_;		
		LineDensityProfile* LineDensityProfile_; 
};
#endif
