/************************************************************************* 
 *
 *    bksolve.c -  all the routines used in the actual numerical solution
 *
 *************************************************************************
 *
 *    This file is part of BKsolver.
 *    Copyright (C) 2004-2005 by Rikard Enberg <REnberg@lbl.gov>
 * 
 *    BKsolver is free software; you can redistribute it and/or 
 *    modify it under the terms of the GNU General Public License 
 *    as published by the Free Software Foundation; either version 2
 *    of the License, or (at your option) any later version.
 * 
 *    BKsolver is distributed in the hope that it will be useful, 
 *    but WITHOUT ANY WARRANTY; without even the implied warranty 
 *    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
 *    See the GNU General Public License for more details.
 * 
 *    You should have received a copy of the GNU General Public 
 *    License along with BKsolver; if not, write to the 
 *    Free Software Foundation, Inc., 59 Temple Place, Suite 330, 
 *    Boston, MA  02111-1307  USA
 *
 *************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_sf_exp.h>
#include "bk.h"

double cheb_int      (int);
double integrand     (int, int);
double cheb_approx   (int);
int    func          (double, const double y[], double f[], void *);
void   evolve        (void);
double hchebint      (int n);
void   realft        (double data[], unsigned long n, int isign);
void   dfct          (int n, double *a, double *t, int *ip, double *w);
void   fftinit       (void);
void   print_grid    (void);
void   print_split   (void);
void   grid_free     (void);

double *twork;
double *wwork;
int    *ipwork;


/*********************************************************************
 *    fftinit: Initialize FFT workspace
 *********************************************************************/
void
fftinit(void)
{
   int iplen = 2+(int)sqrt((chMax+1)/4.); 
   if(((twork = (double *) malloc((chMax+1) * sizeof(double))) == NULL) || 
      ((wwork = (double *) malloc((chMax+1) * sizeof(double))) == NULL) ||
      ((ipwork =   (int *) malloc(iplen * sizeof(int))) == NULL)) 
   {
      fprintf(stderr, "Error: malloc failed in fftinit()");
      exit(EXIT_FAILURE);   
   }
   
   /* This needs to be done first time only, to set up cos/sin table */
   ipwork[0] = 0; 
   
   return;
}

/*********************************************************************
 *    integrand: Return the complete integrand of the integral kernel.
 *               The nonlinear term is defined in cheb_approx.
 *********************************************************************/
double
integrand(int vi, int ui)
{
   double k1, k2, t1, temp, ccfactor;
   double upcut = lnKm2, locut = lnKm1;

   switch (bkKernel) {
      /* Standard BFKL */
      case 0:                 
         t1 = gsl_sf_exp((upcut+locut)*(ugrid[ui]-ugrid[vi]));

         if (ui==vi) {
            k1 = 0.0;
         }
         else {
            k1 = 1.0 / fabs(t1-1.0);
         }
         k2 = 1.0 / sqrt(1 + 4.0/(t1*t1));

         temp  = k1 * (ws[vi] - t1*ws[ui]) + k2 * ws[ui];
         break;
   
      
      /* Two-pole model */
      case 1:  
         if (vi<ui) {
            t1 = gsl_sf_exp((upcut+locut)*(ugrid[vi]-ugrid[ui]));
         }
         else {
            t1=1.0;
         }
         temp = t1 * ws[vi];           
         break;
         
         
      /* BFKL with approximate consistency constraint */
/*       case 2:                 
         t1 = gsl_sf_exp((upcut+locut)*(ugrid[ui]-ugrid[vi]));

         if (vi>ui) {
            ccfactor = exp(0.2*(lnKm2+lnKm1)*(ugrid[ui]-ugrid[vi]));
         } 
         else {
            ccfactor = 1.0;
         }

         if (ui==vi) {
            k1 = 0.0;
         }
         else {
            k1 = 1.0 / fabs(t1-1.0);
         }
         k2 = 1.0 / sqrt(1 + 4.0/(t1*t1));

         temp  = k1 * (ccfactor*ws[vi] - t1*ws[ui]) + k2 * ws[ui];
         break;
 */         
      default:
         temp=0;

   }
   return temp;   
}    



/********************************************************************* 
 *    cheb_int:  Return the value of the integral by approximating 
 *    integrand as Chebyshev series to order chMax. 
 *
 *    Parameters:
 *       u:  grid index of "unintegrated" argument of the integrand
 *
 *    The Chebyshev series for the integrand is computed at 
 *    grid points given by the extrema of the Chebyshev polyomial 
 *    of order chMax. The Chebyshev series for the integral 
 *    is then computed using a fast fourier cosine transform.
 *********************************************************************/
double
cheb_int(int u)
{
   int k;
   double * f;
    
   if((f = (double *) malloc((chMax+1) * sizeof(double))) == NULL) {
      fprintf(stderr, "Error: malloc failed in cheb_int()");
      exit(EXIT_FAILURE);   
   }

   
   /* Fill half the array with the integrand,
      since we only use half the interval - could probably be done better */
   for( k=0 ; k<=chMax/2 ; k++ ) {
      f[k] = integrand(uvMax-k,u);
      f[chMax-k] = 0.0;
   }
   
   /* Call the FFT cosine transform */
   dfct(chMax, f, twork, ipwork, wwork);
   
   double sum = 0.5 * f[0] * hchebint(0);
   for( k=1 ; k<chMax ; k++ ) {
      sum += f[k] * hchebint(k);
   }
   
   sum *= 2.0 /((double)chMax);
   
   free(f);
   return sum;
}

/********************************************************************* 
 *    hchebint(n): return \int _0^1 T_n(x) dx for use in cheb_int
 *********************************************************************/
double
hchebint(int n)
{
   double result;
   if( n==1 ) {
      result = 0.5;
   } else if( GSL_IS_ODD(n) == 1 ) {
      result = (-1.0 + n*sin(n*M_PI_2))/(-1.0 + n*n);
   } else {
      result = 1.0 / (1.0-n*n);
   }
   return result;
}
   
   


/********************************************************************* 
 *    cheb_approx: Compute the rhs as function of ui,
 *                 by integrating over vi and including 
 *                 the non-linear term.  
 *********************************************************************/
double
cheb_approx(int ui)
{
   double M = lnKm1+lnKm2;
   
   double sum = cheb_int(ui);
      
   sum *= M;
   sum -= gsl_pow_2(ws[ui]);

   if(bkKernel==1) {
      sum += (4.0*(M_LN2-1.0)) * ws[ui];
   }
   
   return sum;
}

     
/********************************************************************* 
 *   func: defines the equation for use in the ODE solver
 *********************************************************************/
int
func(double t, const double y[], double f[], void *params)
{
   int i;
   double dummy1; 
   void * dummy2; 
   dummy1 =t;
   dummy2 =params;
   
   for( i=0 ; i<=uvMax ; i++) {
      ws[i] = y[i];
   }
   for( i=0 ; i<=uvMax ; i++) {
      f[i] = abar[i] * cheb_approx(i);
   }
   return GSL_SUCCESS;
}



/********************************************************************* 
 *   evolve: set up ODE system and solve it.
 *           Use Runge-Kutta because the Jacobian is very difficult
 *           to compute.
 *********************************************************************/
void
evolve (void)
{

   unsigned int uvlen = (unsigned int)(uvMax+1);
   const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45;

   gsl_odeiv_step * s    = gsl_odeiv_step_alloc (T, uvlen);
   gsl_odeiv_control * c = gsl_odeiv_control_y_new (odeAccuracy, 0.0);
   gsl_odeiv_evolve * e  = gsl_odeiv_evolve_alloc (uvlen);

   double alpha = 0.0;
   gsl_odeiv_system sys = {func, NULL, uvlen, &alpha};

   double Y = 0.0;
   double h = odeStepsize; 
     
   double * phi;
   if ((phi = (double *) malloc(uvlen * sizeof(double))) == NULL) {
      fprintf(stderr, "Error: malloc failed in evolve()");
      exit(EXIT_FAILURE);   
   }
      
   int i;
   for (i = 0 ; i<=uvMax ; i++) {
      phi[i] = phigrid[0][i];
   }

   fftinit();
   
   if( printSplit == 1 ) {
      print_split();
   }
   else {
      print_grid();
   }
   
   double Ystep = deltaY;
   for (rapidity = 1 ; rapidity <= yMax ; rapidity++) {
       
      double Yi = rapidity * Ystep;
      while (Y < Yi) {
         int status = gsl_odeiv_evolve_apply (e, c, s, &sys, 
                                              &Y, Yi,
                                              &h, phi);
         if (status != GSL_SUCCESS) {
            fprintf(stderr,"Error in gsl_odeiv_evolve_apply !\n");
            fprintf(stderr,"   GSL errorcode: %d\n",status);
            exit(status);
         }
      }
      rapgrid[rapidity] = Y;         
      
      for (i = 0 ; i<=uvMax ; i++) {
         phigrid[rapidity][i] = phi[i];
      }
      
      if( printSplit == 1 ) {
         print_split();
      }
      else {
         print_grid();
      }
   }

   free(phi);
   free(ipwork);
   free(wwork);
   free(twork); 
   
   grid_free();   

   gsl_odeiv_evolve_free (e);
   gsl_odeiv_control_free (c);
   gsl_odeiv_step_free (s);
    
   return;
}
