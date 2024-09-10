/////////////////////////////////////////////////////////////////////////////
//
// FILE NAME
//   ErrorFunctions.cc
//
//   02/19/2015
//
// DESCRIPTION
//   Error functions
//
/////////////////////////////////////////////////////////////////////////////

#include "BufferStore.hh"
#include "OrbitConst.hh"
#include "CppPyWrapper.hh"
#include "ErrorFunctions.hh"

#include <iostream>
#include <cmath>
#include <cfloat>
#include <complex>
#include <typeinfo>

#include "CppPyWrapper.hh"

#include "BufferStore.hh"
#include "OrbitConst.hh"

#include <cmath>
#include <cfloat>
#include <complex>
#include <cstdlib>

// #include <typeinfo>

namespace OrbitUtils{

	std::complex<double> cerrf(std::complex<double> z)
	{
		/** 
		this function calculates the double precision complex error function based on the 
		algorithm of the FORTRAN function written at CERN by K. Koelbig, Program C335, 1970.
	 
		See also M. Bassetti and G.A. Erskine, "Closed expression for the electric field of a 
		two-dimensional Gaussian charge density", CERN-ISR-TH/80-06;
		*/
	
		int n, nc, nu;
		double constant   = 1.12837916709551;
		double xLim = 5.33;
		double yLim = 4.29;
		double h, q, Saux, Sx, Sy, Tn, Tx, Ty, Wx, Wy, xh, xl, x, yh, y;
		double Rx [33];
		double Ry [33];

		x = abs(real(z));
		y = abs(imag(z));
	
		if (y < yLim && x < xLim){
			q = (1.0 - y / yLim) * sqrt(1.0 - (x / xLim) * (x / xLim));
			h  = 1.0 / (3.2 * q);
			nc = 7 + int(23.0 * q);
			xl = pow(h, 1 - nc);
			xh = y + 0.5 / h;
			yh = x;
			nu = 10 + int(21.0 * q);
			Rx[nu] = 0.;
			Ry[nu] = 0.;
			for (n = nu; n > 0; n--){
				Tx = xh + n * Rx[n];
				Ty = yh - n * Ry[n];
				Tn = Tx*Tx + Ty*Ty;
				Rx[n-1] = 0.5 * Tx / Tn;
				Ry[n-1] = 0.5 * Ty / Tn;
 			}
			Sx = 0.;
			Sy = 0.;
			for (n = nc; n>0; n--){
				Saux = Sx + xl;
				Sx = Rx[n-1] * Saux - Ry[n-1] * Sy;
				Sy = Rx[n-1] * Sy + Ry[n-1] * Saux;
				xl = h * xl;
			};
			Wx = constant * Sx;
			Wy = constant * Sy;
		}
		else{
			xh = y;
			yh = x;
			Rx[0] = 0.;
			Ry[0] = 0.;
			for (n = 9; n>0; n--){
				Tx = xh + n * Rx[0];
				Ty = yh - n * Ry[0];
				Tn = Tx * Tx + Ty * Ty;
				Rx[0] = 0.5 * Tx / Tn;
				Ry[0] = 0.5 * Ty / Tn;
			};
			Wx = constant * Rx[0];
			Wy = constant * Ry[0];
		}
		if (y == 0.) {Wx = exp(-x * x);}
		if (imag(z) < 0.){
			Wx =   2.0 * exp(y * y - x * x) * cos(2.0 * x * y) - Wx;
			Wy = - 2.0 * exp(y * y - x * x) * sin(2.0 * x * y) - Wy;
			if (real(z) > 0.) {Wy = -Wy;}
		}
		else if (real(z) < 0.) {Wy = -Wy;}
		
		//std::cout<<" debug Wx="<<Wx<<", Wy="<<Wy<<std::endl;
		return std::complex<double> (Wx, Wy);
	}
}
