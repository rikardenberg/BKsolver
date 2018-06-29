/************************************************************************* 
 *
 *    bksetup.c - various setup routines
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
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_spline.h>
#include "bk.h"


void   set_initcond  (double, int);
void   print_grid    (void);
void   print_split   (void);
void   grid_alloc    (void);
void   grid_free     (void);
void   load_ic       (char *, double tab1[], double tab2[], unsigned int *);
int    check_power   (int num);
char   *xstrdup      (const char *str);



/* Definitions of the global "grids" for storing rapidity, momentum, 
   strong coupling, and the function phi. Furthermore ws is a workspace array. 
   The allocation/freeing is done in grid_alloc() and grid_free() */
double **phigrid;
double *rapgrid;
double *ugrid;
double *ws;
double *abar;



/*********************************************************************
 *    Read initial condition. Adapted from Stephane Munier.
 *********************************************************************/
void
load_ic(char *file, double tab1[], double tab2[], unsigned int *number)
{
   FILE *stream;

   if ( (stream = fopen(file,"r"))==NULL) {
	   fprintf(stderr,"Error opening initial-condition file %s.\n",file);
	   exit(EXIT_FAILURE);
   }
   fprintf(stdout,"\nReading file %s\n",file);
   *number=0;
   while ( !feof(stream) ) {
      fscanf(stream,"%le  %le \n",&tab1[*number],&tab2[*number]);
      (*number)++;
   }
   fflush(stream);  
   fclose(stream);
   fprintf(stdout,"%d points read\n",*number);
   return;
}



/*********************************************************************
 *    Allocate the grids
 *********************************************************************/
void
grid_alloc(void)
{

   int i;
   phigrid = (double **) malloc((yMax+1) * sizeof(double *));
	if( phigrid == NULL ) {
      fprintf(stderr, "Error: malloc failed in grid_alloc()");
      exit(EXIT_FAILURE);   
   }
       
   for(i = 0 ; i < yMax+1 ; i++) {
		phigrid[i] = (double *) malloc((uvMax+1) * sizeof(double));   
	   if( phigrid == NULL ) {
         fprintf(stderr, "Error: malloc failed in grid_alloc()");
         exit(EXIT_FAILURE);   
      }
   }
   
   rapgrid = (double *) malloc((yMax+1)  * sizeof(double));
   ugrid   = (double *) malloc((uvMax+1) * sizeof(double));
   ws      = (double *) malloc((uvMax+1) * sizeof(double));
   abar    = (double *) malloc((uvMax+1) * sizeof(double));
   if((rapgrid == NULL) || (ugrid == NULL) || (ws == NULL) || (abar == NULL)) {
      fprintf(stderr, "Error: malloc failed in grid_alloc()");
      exit(EXIT_FAILURE);   
   }
   
   return;
}


/*********************************************************************
 *    Free the grids
 *********************************************************************/
void
grid_free(void)
{

   int i;
   for(i = 0 ; i < yMax+1 ; i++)
		free(phigrid[i]);   
	
   free(phigrid);   
   free(rapgrid);  
   free(ugrid);    
   free(ws);       
   free(abar);    
   
   return;
}


/*********************************************************************
 *    set_initcond: allocate and set up ugrid, rapgrid, phigrid, ws,
 *                  the running or fixed coupling abar,
 *                  and provide initial condition for phigrid at Y=Y0. 
 *********************************************************************/
void
set_initcond(double satscale, int icchoice)
{
   grid_alloc();
   rapidity = 0;
   rapgrid[rapidity] = 0.0;
   ugrid[0] = 0.0;
   ugrid[uvMax] = 1.0;
   phigrid[0][0] = 0.0;
   
   int i;
   double k2;

   /* Read initial conditions from file into a table and interpolate */
   if (icchoice == 0) {
      unsigned int arrsize;
      double kval;

      double * kt;
      double * ic;
   
      if ((kt = (double *) malloc((10000) * sizeof(double))) == NULL) {
         fprintf(stderr, "Error: malloc failed in set_initcond()");
         exit(EXIT_FAILURE);   
      }
      if ((ic = (double *) malloc((10000) * sizeof(double))) == NULL) {
         fprintf(stderr, "Error: malloc failed in set_initcond()");
         exit(EXIT_FAILURE);   
      }

      load_ic(icFile, kt, ic, &arrsize);

      gsl_interp_accel *acc = gsl_interp_accel_alloc ();
      gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, arrsize);
      gsl_spline_init (spline, kt, ic, arrsize);

      for(i=0 ; i<=uvMax ; i++) {
         ugrid[i] = cos((double)(uvMax-i) * M_PI / chMax);
         kval = gsl_sf_exp(0.5*((lnKm1+lnKm2)*ugrid[i]-lnKm1));
         phigrid[0][i] = gsl_spline_eval (spline, kval, acc);
         
         if(phigrid[0][i]<0.0) phigrid[0][i]=0.0;
      }   
       
      free(kt);
      free(ic); 
      gsl_spline_free (spline);
      gsl_interp_accel_free (acc);
      
   }
   /* Otherwise use simple initial condition */
   else {
      for(i=0 ; i<=uvMax ; i++) {
         ugrid[i] = cos((double)(uvMax-i) * M_PI / chMax);
         k2 = exp((lnKm1+lnKm2)*ugrid[i]-lnKm1);

         switch (icchoice) {
            /* Step function at saturation scale */
            case 1:                 
               if (k2 < satscale)
                  phigrid[0][i] = 1.0;
               else
                  phigrid[0][i] = 0.0;
               break;

            /* Smoothened step function inspired by GBW */
            case 2:          
               if (k2 < satscale)
                  phigrid[0][i] = 1.0-gsl_sf_exp(-gsl_pow_2(k2 - satscale)/4.0);
               else
                  phigrid[0][i] = 0.0;
               break;   

            /* Gaussian in k*phi at saturation scale */
            case 3:     
               gsl_set_error_handler_off();               
               if(i!=uvMax) {
                  phigrid[0][i] = 1.0/sqrt(k2) 
                     * gsl_sf_exp(-gsl_pow_2(log(k2/(satscale*satscale))));
               } 
               else
                  phigrid[0][i] = 0.0;
               gsl_set_error_handler(NULL);
               break;

            /* Not-steep-enough function */
            case 4:          
/*              if (k2 < satscale)
                  phigrid[0][i] = 1.0;
               else
 */                  phigrid[0][i] = 1.0/sqrt(sqrt(k2));
               break;   

            default:
               fprintf(stderr,"Invalid choice of initial conditions!\n");
               exit(EXIT_FAILURE);
               break;
         }
      }
   }
   
   /* Running or fixed coupling */
   double nf = nFlavour;
   double loglambda = log(lambdaQcd);
   double prefac = 12.0/(11.0-2.0*nf/3.0);

   for (i=0 ; i<=uvMax ; i++) {
      if (alphaRunning == 0) {
         abar[i] = alphaBar;
      }
      else {
         abar[i] = prefac / ((lnKm1+lnKm2)*ugrid[i]-lnKm1-2.0*loglambda);
         if (abar[i]>alphaFreeze || abar[i] < 0.0) abar[i] = alphaFreeze;
      }
   }        

   fprintf(stdout,"Initial conditions OK.\n");
   fprintf(stdout,"Creating grid with %d points.\n",uvMax);

   return;
}


/*********************************************************************
   print_grid: print results, all to one file
 *********************************************************************/
void 
print_grid (void)
{
   static int file_exists;
   FILE *fp;
   int i;
   char fil[40];
   
   strncpy(fil,outputFile,40);
   strncat(fil,".dat",4);
   
   if (file_exists) {
      if((fp = fopen (fil, "a+")) == NULL) {
	      fprintf(stderr,"Error opening file %s when appending data.\n",fil);
	      exit(EXIT_FAILURE);
      }
   } 
   else {
      if((fp = fopen (fil, "w")) == NULL) {
	      fprintf(stderr,"Error creating file %s.\n",fil);
	      exit(EXIT_FAILURE);
      }
      fprintf (fp,
         "# k               u                 rapidity          phi\n");
      file_exists = 1;
   }
   
   for (i = 0 ; i<=uvMax ; i++) {
      fprintf (fp,"%+3.6e     %+3.6e     %+3.6e     %+3.6e\n", 
      exp(0.5*((lnKm1+lnKm2)*ugrid[i]-lnKm1)), ugrid[i], rapgrid[rapidity], 
      phigrid[rapidity][i]);
   }
   fprintf (fp, "\n");
   fclose(fp);
}


/*********************************************************************
   print_split: print results, separate file for each rapidity
 *********************************************************************/
void 
print_split (void)
{
   FILE *fp;
   int i;
   char fil[40];
   char suffix[40];
   
   sprintf(suffix,"%04.1f",rapgrid[rapidity]);
   strncpy(fil,outputFile,40);
   strncat(fil,"_",1);
   strncat(fil,suffix,6);
   strncat(fil,".dat",4);
   
   if((fp = fopen (fil, "w")) == NULL) {
	   fprintf(stderr,"Error creating file %s.\n",fil);
	   exit(EXIT_FAILURE);
   }
   
   for (i = 0 ; i<=uvMax ; i++) {
      fprintf (fp,"%+3.6e     %+3.6e\n", 
      exp(0.5*((lnKm1+lnKm2)*ugrid[i]-lnKm1)), 
      phigrid[rapidity][i]);
   }
   
   fprintf (fp, "\n");
   fclose(fp);
}


/*********************************************************************
 *    check_power: Make sure that num == 2^(some integer)
 *********************************************************************/
int 
check_power(int num)
{   
   int i,retval;
   
   if(num<2)
      retval = -1;
      
   for(i = 2 ; i < num ; i = i*2);
   
   if( i == num )
      retval = 0;
   else
      retval = -1;
         
   return retval;  

}



/*********************************************************************
 *    xstrdup: clone of strdup which is not kosher. Nor ANSI C.
 *********************************************************************/
char *
xstrdup (const char *str)
{
  register size_t len = strlen (str) + 1;
  register char *new_str = malloc (len);

  if (new_str != NULL)
    return memcpy (new_str, str, len);

  return NULL;
}
