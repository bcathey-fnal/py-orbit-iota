/////////////////////////////////////////////////////////////////////////////
//
// FILE NAME
//   teapotbase.cc
//
// AUTHOR
//   Jeff Holmes, ORNL, jzh@ornl.gov
//   Joshua Abrams, Knox College, jabrams@knox.edu
//   Steven Bunch, University of Tennessee, sbunch2@utk.edu
//
// Modified by Andrei Shishlo
//   12/30/05
//
// Checked by Jeff Holmes
//   02/2012
//
// DESCRIPTION
//   Define elementary functions for different elements
//
/////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
//
// Include files
//
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
//
// Local Functions:
//
///////////////////////////////////////////////////////////////////////////

#include "teapotbase.hh"
#include "OrbitConst.hh"
#include "Bunch.hh"
#include "SyncPart.hh"

#include <complex>
#include <cmath>

namespace teapot_base
{
    static double* factorial = NULL;

    void init_factorial()
    {
        if(factorial == NULL)
        {
            int n = 50;
            factorial = new double[n];
            factorial[0] = 1.0;
            for(int i = 1; i < n; i++)
            {
                factorial[i] = i * factorial[i - 1];
            }
        }
    }

    void delete_factorial()
    {
        delete [] factorial;
    }

///////////////////////////////////////////////////////////////////////////
// NAME
//   rotatexy
//
// DESCRIPTION
//   Rotates particle coordinates
//
// PARAMETERS
//   bunch = reference to the macro-particle bunch
//   anglexy = rotation angle
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void rotatexy(Bunch* bunch, double anglexy)
{
    double xtemp, pxtemp, ytemp, pytemp;
    double cs = cos(anglexy);
    double sn = sin(anglexy);

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        xtemp  = arr[i][0];
        pxtemp = arr[i][1];
        ytemp  = arr[i][2];
        pytemp = arr[i][3];

        arr[i][0] =  cs * xtemp  - sn * ytemp;
        arr[i][1] =  cs * pxtemp - sn * pytemp;
        arr[i][2] =  sn * xtemp  + cs * ytemp;
        arr[i][3] =  sn * pxtemp + cs * pytemp;
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   drifti
//
// DESCRIPTION
//   Drifts a single particle. Length < 0 is allowed.
//
// PARAMETERS
//   bunch = reference to the macro-particle bunch
//   i = particle index
//   length = length of the drift
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void drifti(Bunch* bunch, int i, double length)
{
    double KNL, phifac, dp_p;

    SyncPart* syncPart = bunch->getSyncPart();

    double gamma2i = 1.0 / (syncPart->getGamma() * syncPart->getGamma());
    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    dp_p = arr[i][5] * dp_p_coeff;
    KNL  = 1.0 / (1.0 + dp_p);
    arr[i][0] += KNL * length * arr[i][1];
    arr[i][2] += KNL * length * arr[i][3];
    phifac = (arr[i][1] * arr[i][1] + arr[i][3] * arr[i][3] +
              dp_p * dp_p * gamma2i) / 2.0;
    phifac = (phifac * KNL - dp_p * gamma2i) * KNL;
    arr[i][4] -= length * phifac;
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   drift
//
// DESCRIPTION
//   Drifts a particle bunch. Length < 0 is allowed.
//
// PARAMETERS
//   bunch = reference to the macro-particle bunch
//   length = length of the drift
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void drift(Bunch* bunch, double length)
{
    double KNL, phifac, dp_p;

    SyncPart* syncPart = bunch->getSyncPart();
    
    double v = OrbitConst::c * syncPart->getBeta();
    syncPart->setTime(syncPart->getTime() + length / v);

    double gamma2i = 1.0 / (syncPart->getGamma() * syncPart->getGamma());
    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        dp_p = arr[i][5] * dp_p_coeff;
        KNL  = 1.0 / (1.0 + dp_p);
        arr[i][0] += KNL * length * arr[i][1];
        arr[i][2] += KNL * length * arr[i][3];
        phifac = (arr[i][1] * arr[i][1] + arr[i][3] * arr[i][3] +
                  dp_p * dp_p * gamma2i) / 2.0;
        phifac = (phifac * KNL - dp_p * gamma2i) * KNL;
        arr[i][4] -= length * phifac;
    }
}
	
///////////////////////////////////////////////////////////////////////////
// NAME
//   wrapbunch
//
// DESCRIPTION
//  wraps the particles longitudinally for the case of a ring beam
//
// PARAMETERS
//  bunch = reference to the macro-particle bunch
//	length = length of the ring
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////
	
void wrapbunch(Bunch* bunch, double length)
{
	//coordinate array [part. index][x,xp,y,yp,z,dE]
	double** arr = bunch->coordArr();
	
	for(int i = 0; i < bunch->getSize(); i++)
		{
			if(arr[i][4] < -length/2.0) arr[i][4] += length;
			if(arr[i][4] > length/2.0) arr[i][4] -= length;
		}
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   kick
//
// DESCRIPTION
//   Kicks a particle bunch
//
// PARAMETERS
//   bunch = reference to the macro-particle bunch
//   kx = strength of the horizontal kick in rad
//   ky = strength of the vertical kick in rad
//   kE = strength of the energy kick in GeV
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void kick(Bunch* bunch, double kx, double ky, double kE, int useCharge)
{
    double charge = +1.0;
    if(useCharge == 1) charge = bunch->getCharge();
    double kxc = kx * charge;
    double kyc = ky * charge;
    double kEc = kE * charge;
    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();
    if(kxc != 0.)
    {
        for(int i = 0; i < bunch->getSize(); i++)
        {
            arr[i][1] += kxc;
        }
    }
    if(kyc != 0.)
    {
        for(int i = 0; i < bunch->getSize(); i++)
        {
            arr[i][3] += kyc;
        }
    }
    if(kEc != 0.)
    {
        for(int i = 0; i < bunch->getSize(); i++)
        {
            arr[i][5] += kEc;
        }
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   multpi
//
// DESCRIPTION
//   Gives particle a multipole momentum kick
//
// PARAMETERS
//   bunch =  reference to the macro-particle bunch
//   i = particle index
//   pole = multipole number
//   kl = integrated strength of the kick [m^(-pole)]
//   skew = 0 - normal, 1 - skew
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void multpi(Bunch* bunch, int i, int pole, double kl, int skew, int useCharge)
{
    double charge = +1.0;
    if(useCharge == 1) charge = bunch->getCharge();
    double klc = kl * charge;
    std::complex<double> z, zn;
    double kl1;

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    kl1 = klc / factorial[pole];
    z = std::complex<double>(arr[i][0], arr[i][2]);

    // take power of z to the n
    zn = std::complex<double>(1.0, 0.0);
    for (int k = 0; k < pole; k++)
    {
        zn *= z;
    }

    // MAD Conventions on signs of multipole terms
    if(skew)
    {
        arr[i][1] += kl1 * std::imag(zn);
        arr[i][3] += kl1 * std::real(zn);
    }
    else
    {
        arr[i][1] -= kl1 * std::real(zn);
        arr[i][3] += kl1 * std::imag(zn);
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   multp
//
// DESCRIPTION
//   Gives particles multipole momentum kicks
//
// PARAMETERS
//   bunch =  reference to the macro-particle bunch
//   pole = multipole number 
//   pole = 0 for dipole, pole = 1 for quad, pole = 2 for sextupole, pole = 3 for octupole
//   kl = integrated strength of the kick [m^(-pole)]
//   skew = 0 - normal, 1 - skew
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void multp(Bunch* bunch, int pole, double kl, int skew, int useCharge)
{
    double charge = +1.0;
    if(useCharge == 1) charge = bunch->getCharge();
    double klc = kl * charge;
    std::complex<double> z, zn;
    double kl1;

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    kl1 = klc / factorial[pole];
    
    for(int i = 0; i < bunch->getSize(); i++)
    {
        z = std::complex<double>(arr[i][0], arr[i][2]);

        // take power of z to the n
        zn = std::complex<double>(1.0, 0.0);
        for (int k = 0; k < pole; k++)
        {
            zn *= z;
        }

        // MAD Conventions on signs of multipole terms
        if(skew)
        {
            arr[i][1] += kl1 * std::imag(zn);
            arr[i][3] += kl1 * std::real(zn);
        }
        else
        {
            arr[i][1] -= kl1 * std::real(zn);
            arr[i][3] += kl1 * std::imag(zn);
        }
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   multpfringeIN
//
// DESCRIPTION
//   Hard edge fringe field for a multipole
//
// PARAMETERS
//   bunch  = reference to the macro-particle bunch
//   pole = multipole number
//   kl = multipole strength
//   skew = multipole skew
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void multpfringeIN(Bunch* bunch, int pole, double kl, int skew, int useCharge)
{
    double charge = +1.0;
    if(useCharge == 1) charge = bunch->getCharge();
    double klc = kl * charge;
    std::complex<double> rootm1 = std::complex<double>(0.0, 1.0);

    SyncPart* syncPart = bunch->getSyncPart();

    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());

    int lm1 = pole;
    int l   = pole + 1;
    int lp1 = pole + 2;
    int lp2 = pole + 3;
    std::complex<double> cxlm1 = std::complex<double>(lm1, 0.0);
    std::complex<double> cxl   = std::complex<double>(l  , 0.0);
    std::complex<double> cxlp1 = std::complex<double>(lp1, 0.0);
    std::complex<double> cxlp2 = std::complex<double>(lp2, 0.0);

    double klfactlp1 = klc / (4.0 * factorial[lp1]);

    // MAD Conventions on signs of multipole terms

    std::complex<double> kterm;
    if(skew)
    {
        kterm = std::complex<double>(0.0, klfactlp1);
    }
    else
    {
        kterm = std::complex<double>(klfactlp1, 0.0);
    }

    double dp_p, KNL;

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        double x = arr[i][0];
        double y = arr[i][2];
        std::complex<double> z = std::complex<double>(x, y);
        double px = arr[i][1];
        double py = arr[i][3];
        dp_p = arr[i][5] * dp_p_coeff;
        KNL  = 1.0 / (1.0 + dp_p);

        // take power of z to the lm1, l

        std::complex<double> zlm1 = std::complex<double>(1., 0.);
        for (int k = 0; k < lm1; k++)
        {
            zlm1 = zlm1 * z;
        }
        std::complex<double> zl = zlm1 * z;

        std::complex<double> fxterm   = std::complex<double>(l * x, -lp2 * y);
        std::complex<double> dxfxterm = cxl * (fxterm + z);
        std::complex<double> dyfxterm = rootm1 * (cxl * fxterm - cxlp2 * z);
        std::complex<double> fyterm   = std::complex<double>(l * y, lp2 * x);
        std::complex<double> dxfyterm = cxl * fyterm + rootm1 * cxlp2 * z;
        std::complex<double> dyfyterm = cxl * (rootm1 * fyterm + z);

        std::complex<double> fxcx   = -kterm * zl   * fxterm;
        std::complex<double> dxfxcx = -kterm * zlm1 * dxfxterm;
        std::complex<double> dyfxcx = -kterm * zlm1 * dyfxterm;
        std::complex<double> fycx   = -kterm * zl   * fyterm;
        std::complex<double> dxfycx = -kterm * zlm1 * dxfyterm;
        std::complex<double> dyfycx = -kterm * zlm1 * dyfyterm;

        arr[i][0] -= std::real(fxcx) * KNL;
        arr[i][2] -= std::real(fycx) * KNL;

        double M11 = 1.0 - std::real(dxfxcx) * KNL;
        double M12 =     - std::real(dxfycx) * KNL;
        double M21 =     - std::real(dyfxcx) * KNL;
        double M22 = 1.0 - std::real(dyfycx) * KNL;
        double detM = M11 * M22 - M12 * M21;

        double pxnew = ( M22 * px - M12 * py) / detM;
        double pynew = (-M21 * px + M11 * py) / detM;

        arr[i][1] = pxnew;
        arr[i][3] = pynew;

        arr[i][4] -= (pxnew * std::real(fxcx) + pynew * std::real(fycx)) *
                     KNL * KNL;
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   multpfringeOUT
//
// DESCRIPTION
//   Hard edge fringe field for a multipole
//
// PARAMETERS
//   bunch  = reference to the macro-particle bunch
//   pole = multipole number
//   kl = multipole strength
//   skew = multipole skew
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void multpfringeOUT(Bunch* bunch, int pole, double kl, int skew, int useCharge)
{
    double charge = +1.0;
    if(useCharge == 1) charge = bunch->getCharge();
    double klc = kl * charge;
    std::complex<double> rootm1 = std::complex<double>(0.0, 1.0);

    SyncPart* syncPart = bunch->getSyncPart();

    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());

    int lm1 = pole;
    int l   = pole + 1;
    int lp1 = pole + 2;
    int lp2 = pole + 3;
    std::complex<double> cxlm1 = std::complex<double>(lm1, 0.0);
    std::complex<double> cxl   = std::complex<double>(l  , 0.0);
    std::complex<double> cxlp1 = std::complex<double>(lp1, 0.0);
    std::complex<double> cxlp2 = std::complex<double>(lp2, 0.0);

    double klfactlp1 = klc / (4.0 * factorial[lp1]);

    // MAD Conventions on signs of multipole terms

    std::complex<double> kterm;
    if(skew)
    {
        kterm = std::complex<double>(0.0, klfactlp1);
    }
    else
    {
        kterm = std::complex<double>(klfactlp1, 0.0);
    }

    double dp_p, KNL;

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        double x = arr[i][0];
        double y = arr[i][2];
        std::complex<double> z = std::complex<double>(x, y);
        double px = arr[i][1];
        double py = arr[i][3];
        dp_p = arr[i][5] * dp_p_coeff;
        KNL  = 1.0 / (1.0 + dp_p);

        // take power of z to the lm1, l

        std::complex<double> zlm1 = std::complex<double>(1., 0.);
        for (int k = 0; k < lm1; k++)
        {
            zlm1 = zlm1 * z;
        }
        std::complex<double> zl = zlm1 * z;

        std::complex<double> fxterm   = std::complex<double>(l * x, -lp2 * y);
        std::complex<double> dxfxterm = cxl * (fxterm + z);
        std::complex<double> dyfxterm = rootm1 * (cxl * fxterm - cxlp2 * z);
        std::complex<double> fyterm   = std::complex<double>(l * y, lp2 * x);
        std::complex<double> dxfyterm = cxl * fyterm + rootm1 * cxlp2 * z;
        std::complex<double> dyfyterm = cxl * (rootm1 * fyterm + z);

        std::complex<double> fxcx   = kterm * zl   * fxterm;
        std::complex<double> dxfxcx = kterm * zlm1 * dxfxterm;
        std::complex<double> dyfxcx = kterm * zlm1 * dyfxterm;
        std::complex<double> fycx   = kterm * zl   * fyterm;
        std::complex<double> dxfycx = kterm * zlm1 * dxfyterm;
        std::complex<double> dyfycx = kterm * zlm1 * dyfyterm;

        arr[i][0] -= std::real(fxcx) * KNL;
        arr[i][2] -= std::real(fycx) * KNL;

        double M11 = 1.0 - std::real(dxfxcx) * KNL;
        double M12 =     - std::real(dxfycx) * KNL;
        double M21 =     - std::real(dyfxcx) * KNL;
        double M22 = 1.0 - std::real(dyfycx) * KNL;
        double detM = M11 * M22 - M12 * M21;

        double pxnew = ( M22 * px - M12 * py) / detM;
        double pynew = (-M21 * px + M11 * py) / detM;

        arr[i][1] = pxnew;
        arr[i][3] = pynew;

        arr[i][4] -= (pxnew * std::real(fxcx) + pynew * std::real(fycx)) *
                     KNL * KNL;
    }
}

////////////////////////////
// NAME
//   quad1
//
// DESCRIPTION
//   Quadrupole element one: linear transport matrix
//
// PARAMETERS
//   bunch  = reference to the macro-particle bunch
//   length = length of transport
//   kq = quadrupole field strength [m^(-2)]
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void quad1(Bunch* bunch, double length, double kq, int useCharge)
{
    double charge = +1.0;
    if(useCharge == 1) charge = bunch->getCharge();
    double kqc = kq * charge;
    if(kqc == 0.)
    {
        drift(bunch,length);
        return;
    }
    double x_init, xp_init, y_init, yp_init;
    double sqrt_kq, kqlength;
    double cx, sx, cy, sy;
    double m11 = 0., m12 = 0., m21 = 0., m22 = 0.;
    double m33 = 0., m34 = 0., m43 = 0., m44 = 0.;

    SyncPart* syncPart = bunch->getSyncPart();

    double v = OrbitConst::c * syncPart->getBeta();
    if(length > 0.)
    {
        syncPart->setTime(syncPart->getTime() + length / v);
    }

    double gamma2i = 1.0 / (syncPart->getGamma() * syncPart->getGamma());
    double dp_p_coeff = 1.0 /(syncPart->getMomentum() * syncPart->getBeta());

    if(kqc > 0.)
    {
        sqrt_kq  = pow(kqc, 0.5);
        kqlength = sqrt_kq * length;
        cx = cos(kqlength);
        sx = sin(kqlength);
        cy = cosh(kqlength);
        sy = sinh(kqlength);
        m11 = cx;
        m12 = sx / sqrt_kq;
        m21 = -sx * sqrt_kq;
        m22 = cx;
        m33 = cy;
        m34 = sy / sqrt_kq;
        m43 = sy * sqrt_kq;
        m44 = cy;
    }
    else if(kqc < 0.)
    {
        sqrt_kq  = pow(-kqc, 0.5);
        kqlength = sqrt_kq * length;
        cx = cosh(kqlength);
        sx = sinh(kqlength);
        cy = cos(kqlength);
        sy = sin(kqlength);
        m11 = cx;
        m12 = sx / sqrt_kq;
        m21 = sx * sqrt_kq;
        m22 = cx;
        m33 = cy;
        m34 = sy / sqrt_kq;
        m43 = -sy * sqrt_kq;
        m44 = cy;
    }

    double dp_p;

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        dp_p    = arr[i][5] * dp_p_coeff;
        x_init  = arr[i][0];
        xp_init = arr[i][1];
        y_init  = arr[i][2];
        yp_init = arr[i][3];

        arr[i][0]  = x_init * m11 + xp_init * m12;
        arr[i][1]  = x_init * m21 + xp_init * m22;
        arr[i][2]  = y_init * m33 + yp_init * m34;
        arr[i][3]  = y_init * m43 + yp_init * m44;
        arr[i][4] += dp_p * gamma2i * length;
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   quad2
//
// DESCRIPTION
//   Quadrupole element two: nonlinear piece
//
// PARAMETERS
//   bunch  = reference to the macro-particle bunch
//   length = length of the element
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void quad2(Bunch* bunch, double length)
{
    double KNL, phifac;

    SyncPart* syncPart = bunch->getSyncPart();

    double gamma2i = 1.0 / (syncPart->getGamma() * syncPart->getGamma());
    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());
    double dp_p;

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        dp_p = arr[i][5] * dp_p_coeff;
        KNL = 1.0 / (1.0 + dp_p);

        arr[i][0] -= KNL * length * dp_p * arr[i][1];
        arr[i][2] -= KNL * length * dp_p * arr[i][3];
        phifac = (arr[i][1] * arr[i][1] + arr[i][3] * arr[i][3] +
                  dp_p * dp_p * gamma2i) / 2.0;
        phifac = (phifac * KNL + dp_p * dp_p * gamma2i) * KNL;
        arr[i][4] -= length * phifac;
    }
}

////////////////////////////
// NAME
//   quad3
//
// DESCRIPTION
//   Quadrupole element 3: non-linear transport 
//   with the longitudinal field component
//
//  It is empty here in the TEAPOT package!
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void quad3(Bunch* bunch, double length, double kq, int useCharge)
{
	return;
}


///////////////////////////////////////////////////////////////////////////
// NAME
//   quadfringeIN
//
// DESCRIPTION
//   Hard edge fringe field for a quad
//
// PARAMETERS
//   bunch =  reference to the macro-particle bunch
//   kq  = strength of quad
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void quadfringeIN(Bunch* bunch, double kq, int useCharge)
{
    double charge = +1.0;
    if(useCharge == 1) charge = bunch->getCharge();
    double kqc = kq * charge;
    double KNL, x_init, xp_init, y_init, yp_init, detM;

    SyncPart* syncPart = bunch->getSyncPart();

    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());
    double dp_p;

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        dp_p    = arr[i][5] * dp_p_coeff;
        KNL     = 1.0 / (1.0 + dp_p);
        x_init  = arr[i][0];
        xp_init = arr[i][1];
        y_init  = arr[i][2];
        yp_init = arr[i][3];
        detM = 1.0 - pow(((kqc * KNL / 4.) *
                          (x_init * x_init - y_init * y_init)), 2);


        arr[i][0] += (kqc * KNL / 12.) * x_init *
                     (x_init * x_init + 3. * y_init * y_init);

        arr[i][1] -= (kqc * KNL / 4.) *
                     (xp_init * (x_init * x_init + y_init * y_init) -
                      2. * yp_init * x_init * y_init);
        arr[i][1] /= detM;

        arr[i][2] -= (kqc * KNL / 12.) * y_init *
                     (y_init * y_init + 3. * x_init * x_init);

        arr[i][3] -= (kqc * KNL / 4.) *
                     (-yp_init * (x_init * x_init + y_init * y_init) +
                      2. * xp_init * x_init * y_init);
        arr[i][3] /= detM;

        arr[i][4] += (kqc * KNL * KNL / 12.) *
                     (xp_init * x_init *
                      (x_init * x_init + 3. * y_init * y_init) -
                      yp_init * y_init *
                      (y_init * y_init + 3. * x_init * x_init));
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   quadfringeOUT
//
// DESCRIPTION
//   Hard edge fringe field for a quad
//
// PARAMETERS
//   bunch  = reference to the macro-particle bunch
//   kq  = strength of quad
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void quadfringeOUT(Bunch* bunch, double kq, int useCharge)
{
    double charge = +1.0;
    if(useCharge == 1) charge = bunch->getCharge();
    double kqc = kq * charge;
    double KNL, x_init, xp_init, y_init, yp_init, detM;

    SyncPart* syncPart = bunch->getSyncPart();

    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());
    double dp_p;

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        dp_p    = arr[i][5] * dp_p_coeff;
        KNL     = 1.0 / (1.0 + dp_p);
        x_init  = arr[i][0];
        xp_init = arr[i][1];
        y_init  = arr[i][2];
        yp_init = arr[i][3];
        detM = 1.0 - pow(((kqc * KNL / 4.) *
                          (x_init * x_init - y_init * y_init)), 2);

        arr[i][0] -= (kqc * KNL / 12.) * x_init *
                     (x_init * x_init + 3. * y_init * y_init);

        arr[i][1] += (kqc * KNL / 4.) *
                     (xp_init * (x_init * x_init + y_init * y_init) -
                      2. * yp_init * x_init * y_init);
        arr[i][1] /= detM;

        arr[i][2] += (kqc * KNL / 12.) * y_init *
                     (y_init * y_init + 3. * x_init * x_init);

        arr[i][3] += (kqc * KNL / 4.) *
                     (-yp_init * (x_init * x_init + y_init * y_init) +
                      2. * xp_init * x_init * y_init);
        arr[i][3] /= detM;

        arr[i][4] -= (kqc * KNL * KNL / 12.) *
                     (xp_init * x_init *
                      (x_init * x_init + 3. * y_init * y_init) -
                      yp_init * y_init *
                      (y_init * y_init + 3. * x_init * x_init));
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   wedgerotate
//
// DESCRIPTION
//   Rotates coordinates by e for fringe fields at non-SBEND
//
// PARAMETERS
//   bunch  = reference to the macro-particle bunch
//   e = rotation angle
//   frinout = 0 before fringe, 1 after fringe
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void wedgerotate(Bunch* bunch, double e, int frinout)
{
    double cs, sn;
    double xp_temp, p0_temp, p0;

    SyncPart* syncPart = bunch->getSyncPart();

    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());
    double dp_p;

    cs = cos(e);
    sn = sin(e);

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        if(frinout == 0)
        {
            dp_p    = arr[i][5] * dp_p_coeff;
            xp_temp = arr[i][1];
            p0_temp = 1.0 + dp_p;

            arr[i][0] /=  cs;
            arr[i][1]  =  xp_temp * cs + p0_temp * sn;
            p0         = -xp_temp * sn + p0_temp * cs;
            dp_p       =  p0 - 1.0;
            arr[i][4]  =  (-arr[i][0] * sn + arr[i][4]) * cs;
        }
        else
        {
            dp_p = arr[i][5] * dp_p_coeff;
            p0   = 1.0 + dp_p;

            arr[i][4]  = arr[i][0] * sn + arr[i][4] / cs;
            arr[i][0] *= cs;
            xp_temp    = arr[i][1] * cs - p0 * sn;
            p0_temp    = arr[i][1] * sn + p0 * cs;
            arr[i][1]  = xp_temp;
            dp_p       = p0_temp - 1.0;
        }
        arr[i][5] = dp_p / dp_p_coeff;
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   wedgedrift
//
// DESCRIPTION
//   Drifts particles through wedge for non-SBEND
//
// PARAMETERS
//   bunch =  reference to the macro-particle bunch
//   e = wedge angle
//   inout = 0 for in, 1 for out
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void wedgedrift(Bunch* bunch, double e, int inout)
{
    double ct, tn;
    double s;

    SyncPart* syncPart = bunch->getSyncPart();

    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());
    double dp_p;

    ct = cos(e) / sin(e);

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        if(inout == 0)
        {
            dp_p = arr[i][5] * dp_p_coeff;
            tn   = arr[i][1] / (1.0 + dp_p);
            s    = arr[i][0] / (ct - tn);
        }
        else
        {
            s    = arr[i][0] / ct;
        }

        drifti(bunch, i, s);
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   wedgebend
//
// DESCRIPTION
//   Straight bends particles through wedge for non-SBEND
//
// PARAMETERS
//   bunch  = reference to the macro-particle bunch
//   e = wedge angle
//   inout = 0 for in, 1 for out
//   rho = radius of curvature
//   nsteps = number of integraton steps
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void wedgebend(Bunch* bunch, double e, int inout, double rho, int nsteps)
{
    double ct, tn;
    double s, sm, sm2;
    int nst;

    SyncPart* syncPart = bunch->getSyncPart();

    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());
    double dp_p;

    nst = nsteps / 2;
    if(nst < 1) nst = 1;
    ct = cos(e) / sin(e);

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        if(inout == 0)
        {
            s    = -arr[i][0] / ct;
        }
        else
        {
            dp_p =  arr[i][5] * dp_p_coeff;
            tn   =  arr[i][1] / (1.0 + dp_p);
            s    = -arr[i][0] / (ct + tn);
        }

        sm  = s / nst;
        sm2 = sm / 2.0;

        drifti(bunch, i, sm2);
        arr[i][1] -= sm / rho;
        for(int j  = 1; j < nst; j++)
        {
            drifti(bunch, i, sm);
            arr[i][1] -= sm / rho;
        }
        drifti(bunch, i, sm2);
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   bend1
//
// DESCRIPTION
//   Linear bend transport
//
// PARAMETERS
//   bunch  = reference to the macro-particle bunch
//   length = length of transport
//   th = bending angle
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void bend1(Bunch* bunch, double length, double th)
{
    double x_init, xp_init;
    double cx, sx, rho;
    double m11, m12, m16;
    double m21, m22, m26;
    double m51, m52, m56;

    SyncPart* syncPart = bunch->getSyncPart();

    double v = OrbitConst::c * syncPart->getBeta();
    if(length > 0.)
    {
	   syncPart->setTime( syncPart->getTime() + length/v);
    }

    double betasq = syncPart->getBeta() * syncPart->getBeta();
    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());
    double dp_p;

    rho = length / th;
    cx  = cos(th);
    sx  = sin(th);
    m11 = cx;
    m12 = rho * sx;
    m16 = rho * (1.0 - cx);
    m21 = -sx / rho;
    m22 = cx;
    m26 = sx;
    m51 = -sx;
    m52 = -rho * (1.0 - cx);
    m56 = -betasq * length + rho * sx;

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        dp_p   = arr[i][5] * dp_p_coeff;
        x_init = arr[i][0];
        xp_init = arr[i][1];

        arr[i][0]  = x_init * m11 + xp_init * m12 + dp_p * m16;
        arr[i][1]  = x_init * m21 + xp_init * m22 + dp_p * m26;
        arr[i][2] += length * arr[i][3];
        arr[i][4] += x_init * m51 + xp_init * m52 + dp_p * m56;
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   bend2
//
// DESCRIPTION
//   Kinetic bend transport (same as nonlinear quad transport - quad2)
//
// PARAMETERS
//   bunch  = reference to the macro-particle bunch
//   length = length of element (either full of half step)
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void bend2(Bunch* bunch, double length)
{
    double KNL, phifac;

    SyncPart* syncPart = bunch->getSyncPart();

    double gamma2i = 1.0 / (syncPart->getGamma() * syncPart->getGamma());
    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());
    double dp_p;

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        dp_p = arr[i][5] * dp_p_coeff;
        KNL = 1.0 / (1.0 + dp_p);

        arr[i][0] -= KNL * length * dp_p * arr[i][1];
        arr[i][2] -= KNL * length * dp_p * arr[i][3];
        phifac = (arr[i][1] * arr[i][1] + arr[i][3] * arr[i][3] +
                  dp_p * dp_p * gamma2i) / 2.0;
        phifac = (phifac * KNL + dp_p * dp_p * gamma2i) * KNL;
        arr[i][4] -= length * phifac;
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   bend3
//
// DESCRIPTION
//   Nonlinear curvature bend transport
//   depending on py and dE in Hamiltonian
//
// PARAMETERS
//   bunch =  reference to the macro-particle bunch
//   th = bending angle
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void bend3(Bunch* bunch, double th)
{
    double KNL, phifac;

    SyncPart* syncPart = bunch->getSyncPart();

    double gamma2i = 1.0 / (syncPart->getGamma() * syncPart->getGamma());
    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());
    double dp_p;

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        dp_p = arr[i][5] * dp_p_coeff;
        KNL  = 1.0 / (1.0 + dp_p);

        phifac     = (arr[i][3] * arr[i][3] + dp_p * dp_p * gamma2i) / 2.0;
        arr[i][1] -= phifac * KNL * th;
        arr[i][2] += KNL * arr[i][3] * arr[i][0] * th;
        phifac     = (phifac * KNL - dp_p * gamma2i) * KNL;
        arr[i][4] -= th * phifac * arr[i][0];
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   bend4
//
// DESCRIPTION
//   Nonlinear curvature bend transport
//   depending on px in Hamiltonian
//
// PARAMETERS
//   bunch =  reference to the macro-particle bunch
//   th = bending angle
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void bend4(Bunch* bunch, double th)
{
    double KNL, phifac, xfac;

    SyncPart* syncPart = bunch->getSyncPart();

    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());
    double dp_p;

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        dp_p = arr[i][5] * dp_p_coeff;
        KNL  = 1.0 / (1.0 + dp_p);

        xfac   = 1.0 + KNL * arr[i][1] * th / 2.0;
        phifac = KNL * KNL * arr[i][1] * arr[i][1] / 2.0;
        arr[i][0] *= xfac * xfac;
        arr[i][1] /= xfac;
        arr[i][4] -= th * phifac * arr[i][0];
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   driftexact
//
// DESCRIPTION
//   Exact drift transport. Where drift() uses the Hamiltonian of Eq. (8) of
//   J. Holmes, "Single Particle Transport in ORBIT and pyORBIT" (2022), in
//   which the kinematic square root is expanded to second order in the
//   transverse momenta, this routine keeps the square root itself:
//
//     H = dE - sqrt(P^2 - px^2 - py^2),  P^2 = (1 + dE)^2 - dE^2 / gamma^2
//
//   where dE is the energy deviation scaled by beta^2 * E0 (the "dp_p" of the
//   other routines) and P is the exact momentum in units of the design
//   momentum. This is the th -> 0 limit of bendexact().
//
// PARAMETERS
//   bunch  = reference to the macro-particle bunch
//   length = length of transport
//
// RETURNS
//   Nothing
//
// - nilanjan@fnal.gov, 09/07/2026
///////////////////////////////////////////////////////////////////////////

void driftexact(Bunch* bunch, double length)
{
    SyncPart* syncPart = bunch->getSyncPart();

    double v = OrbitConst::c * syncPart->getBeta();
    if(length > 0.)
    {
        syncPart->setTime(syncPart->getTime() + length / v);
    }

    double gamma2i = 1.0 / (syncPart->getGamma() * syncPart->getGamma());
    double beta2 = 1.0 - gamma2i;
    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());
    double dp_p, w, pm1, pz2, pz;

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        dp_p = arr[i][5] * dp_p_coeff;

        // pm1 = P^2 - 1 and w = d(P^2 / 2) / d(dE), both exact and both
        // written so that they carry no cancellation for a small dE
        pm1 = dp_p * (2.0 + beta2 * dp_p);
        w   = 1.0 + beta2 * dp_p;

        pz2 = 1.0 + pm1 - arr[i][1] * arr[i][1] - arr[i][3] * arr[i][3];
        if(pz2 <= 0.0)
        {
            // the particle is not moving forward - it has no image under this map
            bunch->deleteParticleFast(i);
            continue;
        }
        pz = sqrt(pz2);

        arr[i][0] += length * arr[i][1] / pz;
        arr[i][2] += length * arr[i][3] / pz;
        arr[i][4] += length * (pz - w) / pz;
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   bendexact
//
// DESCRIPTION
//   Exact sector bend transport. This is the closed-form map of the full
//   sector bend Hamiltonian, with the kinematic square root kept rather than
//   expanded:
//
//     H = (r^2 - 1) / 2 + dE - r * sqrt(P^2 - px^2 - py^2),   r = 1 + x / rho
//
//   with P^2 = (1 + dE)^2 - dE^2 / gamma^2 the exact momentum in units of the
//   design momentum. It replaces the combination bend1 + bend2 + bend3 + bend4,
//   which together integrate the same Hamiltonian with the square root
//   expanded to second order in px, py and dE / gamma. Because this map is
//   exact, a pure sector dipole is transported exactly in a single step, with
//   no dependence on the number of steps.
//
//   The motion is a circular arc: py, dE and pperp^2 = P^2 - py^2 are
//   conserved. The map is that of ImpactX ExactCFbend's H_1, rewritten in the
//   pyORBIT variables and regrouped so that no term suffers cancellation as
//   the bend angle goes to zero, in which limit it reduces exactly to
//   driftexact().
//
// PARAMETERS
//   bunch  = reference to the macro-particle bunch
//   length = length of transport
//   th = bending angle
//
// RETURNS
//   Nothing
//
// - nilanjan@fnal.gov, 09/07/2026
///////////////////////////////////////////////////////////////////////////

void bendexact(Bunch* bunch, double length, double th)
{
    // nothing is transported over a vanishing length
    if(fabs(length) < OrbitConst::tiny) return;

    // a bend with a vanishing angle is an exact drift
    if(fabs(th) < OrbitConst::tiny)
    {
        driftexact(bunch, length);
        return;
    }

    SyncPart* syncPart = bunch->getSyncPart();

    double v = OrbitConst::c * syncPart->getBeta();
    if(length > 0.)
    {
        syncPart->setTime(syncPart->getTime() + length / v);
    }

    double gamma2i = 1.0 / (syncPart->getGamma() * syncPart->getGamma());
    double beta2 = 1.0 - gamma2i;
    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());

    double rho = length / th;
    double cx = cos(th);
    double sx = sin(th);
    double sh = sin(0.5 * th);

    // rho * (1 - cos(th)) and rho * sin(th), formed without cancellation so
    // that they stay accurate however small the bend angle of a step becomes
    double rho_omc = 2.0 * rho * sh * sh;
    double rho_sin = length * (sx / th);

    // below one radian per step the two arcsines below are close enough that
    // their difference must be taken through the arcsine addition identity
    int smallAngle = (fabs(th) < 1.0) ? 1 : 0;

    double dp_p, w, pm1, x, px, py2, pperp2, pperp;
    double pzi, pzf, pzi2, pzf2, pzim1, pxout;
    double rho_dpx, rho_dpz, rho_num, u, rho_dasin;

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        dp_p = arr[i][5] * dp_p_coeff;

        // pm1 = P^2 - 1 and w = d(P^2 / 2) / d(dE), both exact
        pm1 = dp_p * (2.0 + beta2 * dp_p);
        w   = 1.0 + beta2 * dp_p;

        x   = arr[i][0];
        px  = arr[i][1];
        py2 = arr[i][3] * arr[i][3];

        // pperp^2 = P^2 - py^2 is a constant of the motion, as are py and dE
        pperp2 = 1.0 + pm1 - py2;
        pzi2   = pperp2 - px * px;
        if(pzi2 <= 0.0)
        {
            // the particle is not moving forward - it has no image under this map
            bunch->deleteParticleFast(i);
            continue;
        }
        pperp = sqrt(pperp2);
        pzi   = sqrt(pzi2);

        // pzi - 1 without cancellation
        pzim1 = (pm1 - py2 - px * px) / (pzi + 1.0);

        pxout = px * cx + (pzim1 - x / rho) * sx;
        pzf2  = pperp2 - pxout * pxout;
        if(pzf2 <= 0.0)
        {
            bunch->deleteParticleFast(i);
            continue;
        }
        pzf = sqrt(pzf2);

        // rho * (px - pxout) and rho * (pzf - pzi), both free of cancellation
        rho_dpx = px * rho_omc - rho_sin * pzim1 + x * sx;
        rho_dpz = rho_dpx * (px + pxout) / (pzf + pzi);

        // rho * [asin(px / pperp) - asin(pxout / pperp)], the arc through
        // which the momentum vector turns, less the design bend angle. The
        // identity asin(a) - asin(b) = asin(a * sqrt(1 - b^2) - b * sqrt(1 - a^2))
        // moves the cancellation inside the arcsine, where
        // rho * (px * pzf - pxout * pzi) = px * rho_dpz + pzi * rho_dpx.
        if(smallAngle)
        {
            rho_num = px * rho_dpz + pzi * rho_dpx;
            u = rho_num / (rho * pperp2);
            rho_dasin = (rho_num / pperp2) *
                        ((fabs(u) < 1.0e-4) ? (1.0 + u * u / 6.0) : (asin(u) / u));
        }
        else
        {
            rho_dasin = rho * (asin(px / pperp) - asin(pxout / pperp));
        }

        // x_out = x cos(th) + rho (1 - cos(th)) (pzi - 1)
        //         + rho (pzf - pzi) + rho px sin(th)
        arr[i][0]  = x * cx + rho_omc * pzim1 + rho_dpz + rho_sin * px;
        arr[i][1]  = pxout;
        arr[i][2] += (length + rho_dasin) * arr[i][3];
        // z picks up length - w * rho * theta, with theta the turned arc; the
        // leading parts of the two terms cancel and are removed analytically
        arr[i][4] += -length * beta2 * dp_p - w * rho_dasin;
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   bendfringeIN
//
// DESCRIPTION
//   Hard edge fringe field for a bend
//
// PARAMETERS
//   bunch =  reference to the macro-particle bunch
//   rho = radius of curvature for bending
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void bendfringeIN(Bunch* bunch, double rho)
{
    SyncPart* syncPart = bunch->getSyncPart();

    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());
    double dp_p, KNL;

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        dp_p    = arr[i][5] * dp_p_coeff;
        KNL  = 1.0 / (1.0 + dp_p);

        arr[i][0] += KNL * arr[i][2] * arr[i][2] / (2. * rho);
        arr[i][3] -= KNL * arr[i][1] * arr[i][2] / rho;
        arr[i][4] -= KNL * KNL * arr[i][1] * arr[i][2] * arr[i][2] / (2. * rho);
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   bendfringeOUT
//
// DESCRIPTION
//   Hard edge fringe field for a bend
//
// PARAMETERS
//   bunch =  reference to the macro-particle bunch
//   rho = radius of curvature for bending
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void bendfringeOUT(Bunch* bunch, double rho)
{
    SyncPart* syncPart = bunch->getSyncPart();

    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());
    double dp_p, KNL;

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        dp_p    = arr[i][5] * dp_p_coeff;
        KNL  = 1.0 / (1.0 + dp_p);

        arr[i][0] -= KNL * arr[i][2] * arr[i][2] / (2. * rho);
        arr[i][3] += KNL * arr[i][1] * arr[i][2] / rho;
        arr[i][4] += KNL * KNL * arr[i][1] * arr[i][2] * arr[i][2] / (2. * rho);
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   soln
//
// DESCRIPTION
//   Integration through a solenoid
//
// PARAMETERS
//   bunch  =  reference to the macro-particle bunch
//   length = integration length
//   B      = magnetic field (1/m)
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void soln(Bunch* bunch, double length, double B, int useCharge)
{
    // A solenoid with zero field is just a drift - nilanjan@fnal.gov, 02/04/24
    if(abs(B) < OrbitConst::tiny)
    {
        drift(bunch, length);
        return;
    }

    double charge = +1.0;
    if(useCharge == 1) charge = bunch->getCharge();
    double Bc = B * charge;
    double KNL, phase, cs, sn;
    double cu, cpu, u_init, pu_init, u, pu, phifac;

    SyncPart* syncPart = bunch->getSyncPart();

    double v = OrbitConst::c * syncPart->getBeta();
    if(length > 0.)
    {
	   syncPart->setTime( syncPart->getTime() + length/v);
    }

    double gamma2i = 1.0 / (syncPart->getGamma() * syncPart->getGamma());
    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());
    double dp_p;

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        dp_p = arr[i][5] * dp_p_coeff;
        KNL  = 1.0 / (1.0 + dp_p);

        cu      =  arr[i][2] / 2.     - arr[i][1] / Bc;
        cpu     =  arr[i][0] * Bc / 2. + arr[i][3];
        u_init  =  arr[i][2] / 2.     + arr[i][1] / Bc;
        pu_init = -arr[i][0] * Bc / 2. + arr[i][3];
        phase = KNL * Bc * length;
        cs = cos(phase);
        sn = sin(phase);

        u =   u_init * cs     + pu_init * sn / Bc;
        pu = -u_init * Bc * sn + pu_init * cs;

        arr[i][0] = (-pu + cpu) / Bc;
        arr[i][1] = 0.5 * (u - cu) * Bc;
        arr[i][2] = u + cu;
        arr[i][3] = 0.5 * (pu + cpu);

        phifac = (pu_init * pu_init +
                  Bc * Bc * u_init * u_init +
                  dp_p * dp_p * gamma2i
                 ) / 2.0;
        phifac = (phifac * KNL - dp_p * gamma2i) * KNL;
        arr[i][4] -= length * phifac;
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   wedgebendCF
//
// DESCRIPTION
//   Straight bends particles through wedge for Combined Function non-SBEND
//
// PARAMETERS
//   bunch =  reference to the macro-particle bunch
//   e = wedge angle
//   inout = 0 for in, 1 for out
//   rho = radius of curvature
//   vecnum = number of multipole terms
//   pole = multipolarities of multipole terms
//   kl = integrated strengths of multipole terms
//   skew = skewness  of multipole terms
//   nsteps = number of integraton steps
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void wedgebendCF(Bunch* bunch, double e, int inout,
                 double rho,
                 int vecnum,
                 std::vector<int>& pole,
                 std::vector<double>& kl,
                 std::vector<int>& skew,
                 int nsteps, int useCharge)
{
    double ct, tn;
    double s, sm, sm2, klint;
    int nst;

    SyncPart* syncPart = bunch->getSyncPart();

    double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());
    double dp_p;

    nst = nsteps / 2;
    if(nst < 1) nst = 1;
    ct = cos(e) / sin(e);

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        if(inout == 0)
        {
            s = -arr[i][0] / ct;
        }
        else
        {
            dp_p =  arr[i][5] * dp_p_coeff;
            tn   =  arr[i][1] / (1.0 + dp_p);
            s    = -arr[i][0] / (ct + tn);
        }

        sm = s / nst;
        sm2 = sm / 2.0;

        drifti(bunch, i, sm2);
        arr[i][1] -= sm / rho;
        for (int l = 0; l < vecnum; l++)
        {
            klint = kl[l] * sm;
            multpi(bunch, i, pole[l], klint, skew[l], useCharge);
        }
        for(int j = 1; j < nst; j++)
        {
            drifti(bunch, i, sm);
            arr[i][1] -= sm / rho;
            for (int l = 0; l < vecnum; l++)
            {
                klint = kl[l] * sm;
                multpi(bunch, i, pole[l], klint, skew[l], useCharge);
            }
        }
        drifti(bunch, i, sm2);
    }
}

///////////////////////////////////////////////////////////////////////////
// NAME
//   RingRF
//
// DESCRIPTION
//   Ring type RF cavity. Transition time factor T(k) = const = T(k0).
//   No need for symplectic phase correction.
//
// PARAMETERS
//   bunch =  reference to the macro-particle bunch
//   harmonic_numb = harmonic number
//   voltage = voltage in Giga Volts
//   phase_s = synchronous phase in Rad
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void RingRF(Bunch* bunch, double ring_length, int harmonic_numb,
            double voltage, double phase_s, int useCharge)
{
    double charge = +1.0;
    if(useCharge == 1) charge = bunch->getCharge();
    double deltaV = 0.;
    double coeff  = charge;

    double Factor = 0.;
    if(ring_length > 0.)
    {
        Factor = 2.0 * OrbitConst::PI/ring_length;
    }

    SyncPart* syncPart = bunch->getSyncPart();

    if(phase_s != 0.)
    {
        double kin_e = syncPart->getEnergy();
        kin_e += coeff * voltage * sin(phase_s);
        syncPart->setMomentum(syncPart->energyToMomentum(kin_e));
    }

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        deltaV = voltage * ( sin(harmonic_numb*Factor*arr[i][4] + phase_s));
        arr[i][5] += coeff * deltaV;
    }
}

// New element added - nilanjan@uchicago.edu, 10/18/21
// Based on https://arxiv.org/pdf/2106.03327v2.pdf
// See also SLAC-75, pg 117
///////////////////////////////////////////////////////////////////////////
// NAME
//   dipedge
//
// DESCRIPTION
//   Zero length linear dipole edge transport (same as MAD-X dipedge)
//
// PARAMETERS
//   bunch  = reference to the macro-particle bunch
//   h      = angle/length - must equal the associated sbend to make sense
//   e1     = edge angle
//   fint   = field integral - same definition as for a MAD-X SBEND
//   hgap   = half gap height of the associated sbend
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

void dipedge(Bunch* bunch, double h, double e1, double fint, double hgap)
{
    SyncPart* syncPart = bunch->getSyncPart();

    // double dp_p_coeff = 1.0 / (syncPart->getMomentum() * syncPart->getBeta());
    double sin_e1 = sin(e1);
    double h_tan_e1 = h*tan(e1);
    double psi = 2.0*fint*hgap*h*(1.0+sin_e1*sin_e1)/cos(e1);
    double h_tan_e1_minus_psi = h*tan(e1-psi);

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();

    for(int i = 0; i < bunch->getSize(); i++)
    {
        //double dp_p = arr[i][5] * dp_p_coeff;

        arr[i][1] += h_tan_e1*arr[i][0];
        arr[i][3] -= h_tan_e1_minus_psi*arr[i][2];

    }
}


}  //end of namespace teapot_base
