/********************************************************************* 
 *
 *    Calculate saturation momentum, its Y-derivative, and the 
 *    reduced front shape from a previously generated grid
 *    Usage: qs <filename> 
 *           where <filename> contains the (k,Y)-grid in the 
 *           format output by the bkeqn.c program
 *            
 *********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>

#define UVMAX 2000
#define YMAX 300

void   read_grid  (char *, int *, int *, double *, double *);
double phifunc    (double, int);
double bisect     (double, double, double, int);
double bisect2    (double, double, double, int);
double dfunc      (double, double *);
double func       (double, void *);
double dtest      (double, void *);
double analytic   (double, int);
double reddqs     (double, double, int);
double reddqs2    (double, double, int);

double phigrid[YMAX+1][UVMAX+1];
double ugrid[UVMAX+1];
double kgrid[UVMAX+1];
double rapgrid[YMAX+1];

int    umax,rapmax;
double lnk2min,lnk2max;

int main (int argc, char **argv) 
{

   double qs[YMAX+1];
   double logqs2[YMAX+1];
   int yi;
   int irap;
   double rap;
   double deriv;
   double deltay;
   double gammac = 0.6275;
    
   FILE *fp1, *fp2, *fp3, *fp4;
   char file_redfr[100];   /* Reduced front                         */
   char file_qs[100];      /* Saturation scale                      */
   char file_dqs[100];     /* Derivative of sat. scale              */
   char file_q3rd[100];    /* Test of the third term (reduced dQs)  */
   char suffix[100] = "";  /* Optional filename suffix              */

/*--------------------------------------------------------*/
/*      Take care of setting up filenames etc             */
/*--------------------------------------------------------*/
   
   
   if( argc<2 ) {
      printf("Usage: %s <grid-file> [<optional name suffix>]\n",argv[0]); 
      printf("       where <grid-file> contains the grid in correct format\n");
      printf("       and <optional name suffix> will be added to filename.\n");
      printf("Output files:\n");
      printf("       qs_<grid-file>    for Qs(Y)\n");
      printf("       dqs_<grid-file>   for dlog Qs(Y)^2/dY\n");
      printf("       redfr_<grid-file> for reduced front.\n");
      printf("       q3rd_<grid-file>  for reduced dQs.\n");
      printf("Files will be overwritten if they exist.\n");
      return 1;
   }
   
   if( argc>3 ) {
	   fprintf(stderr,"Warning: extra command line arguments ignored!\n");
   }

   if ( argc>2 ) {
      sprintf(suffix,"%s_",argv[2]);
   } 
     
   
   sprintf(file_redfr,"redfr_%s%s",suffix,argv[1]);
   sprintf(file_qs,"qs_%s%s",suffix,argv[1]);
   sprintf(file_dqs,"dqs_%s%s",suffix,argv[1]);
   sprintf(file_q3rd,"q3rd_%s%s",suffix,argv[1]);
    
   if ( (fp1 = fopen (file_qs, "w")) == NULL ) {
	   fprintf(stderr,"Error creating file %s\n",file_qs);
	   exit(EXIT_FAILURE);
   }
   if ( (fp2 = fopen (file_dqs, "w")) == NULL ) {
	   fprintf(stderr,"Error creating file %s\n",file_dqs);
	   exit(EXIT_FAILURE);
   }
   if ( (fp3 = fopen (file_redfr, "w")) == NULL ) {
      fprintf(stderr,"Error creating file %s\n",file_redfr);
      exit(EXIT_FAILURE);
   }
   if ( (fp4 = fopen (file_q3rd, "w")) == NULL ) {
      fprintf(stderr,"Error creating file %s\n",file_q3rd);
      exit(EXIT_FAILURE);
   }
	printf("Writing Qs to file:            %s\n",file_qs);
	printf("Writing dlogQs2/dY to file:    %s\n",file_dqs);
   printf("Writing reduced front to file: %s\n",file_redfr);
   printf("Writing reduced dQs to file:   %s\n\n",file_q3rd);

   read_grid(argv[1],&rapmax,&umax,&lnk2min,&lnk2max);
   printf("# rapmax = %d, umax = %d, lnk2min = %g, lnk2max = %g\n",
           rapmax,umax,lnk2min,lnk2max);
           
           
           
   
/*--------------------------------------------------------*/
/*      Find the saturation momentum for each rapidity,   */
/*      i.e. where the amplitude = e.g. 0.01              */
/*--------------------------------------------------------*/
   
   for ( yi = 1 ; yi <= rapmax ; yi++ ) {
      qs[yi] = bisect(0.01, 1e-1, 1e60, yi);      
      logqs2[yi] = log(qs[yi]*qs[yi]);      
   }
   
   fprintf(fp1,"# rapidity      qs\n");
   for ( yi = 1 ; yi <= rapmax ; yi++ ) {
      fprintf(fp1,"  %3.6e  %3.6e\n", rapgrid[yi], qs[yi]);
   }
   

  
/*--------------------------------------------------------*/
/*      Compute the derivative dlog(Qs^2)/dY              */
/*      and the reduced sat scale derivative              */
/*--------------------------------------------------------*/
    
   deltay = rapgrid[2]-rapgrid[1];
    
   fprintf(fp2,"# rapidity      num deriv     analytic-3    analytic-2\n");
   fprintf(fp4,"# rapidity      F(Y)          G(Y)\n");
   
   for ( irap = 1 ; irap < rapmax ; irap++ ) {
         
      rap = rapgrid[irap];
          
      if (irap<=2) {
        deriv = (0.5/deltay) *
           (-logqs2[irap+2] + 4.0*logqs2[irap+1] - 3.0*logqs2[irap]);
      } 
      else if (irap>=rapmax-2) {
        deriv = (0.5/deltay) *
           (logqs2[irap-2] - 4.0*logqs2[irap-1] + 3.0*logqs2[irap]);
      } 
      else {           
        deriv = (1.0/(12.0*deltay)) *
           (-logqs2[irap+2] + 8.0*logqs2[irap+1]
            -8.0*logqs2[irap-1] + logqs2[irap-2]);
      } 
          
      fprintf(fp2,"  %3.6e  %3.6e  %3.6e  %3.6e\n", 
        rap, deriv, analytic(rap,3), analytic(rap,2));
           
      fprintf(fp4,"  %3.6e  %3.6e  %3.6e  %3.6e\n",
           rap, reddqs2(rap,deriv,3), reddqs2(rap,deriv,2), 
           reddqs2(rap,analytic(rap,3),2));           
   }



/*--------------------------------------------------------*/
/*      Plot the reduced front                            */
/*--------------------------------------------------------*/
                                                           
   fprintf (fp3,"\n# log[k^2/Qs^2]   rapidity          red. phi\n");
   for ( yi = 1 ; yi <= rapmax ; yi++ ) {
      int i;
      for (i = 0 ; i<=umax ; i++) {
         fprintf (fp3,"%+3.6e     %+3.6e     %+3.6e\n", 
         2.0*log(kgrid[i]/qs[yi]),
         rapgrid[yi],
         exp(gammac*2.0*log(kgrid[i]/qs[yi]))*(phigrid[yi][i])); 
      }
      fprintf(fp3,"\n");
   }
    


/*--------------------------------------------------------*/
/*      Finish                                            */
/*--------------------------------------------------------*/
   fclose(fp1);                                            
   fclose(fp2);
   fclose(fp3);
   fclose(fp4);
   return 0;
}



/*********************************************************************
 *********************************************************************
 *********************************************************************/


/*********************************************************************
 *    func: numerical interpolation of log Qs(Y)^2
 *********************************************************************/
double 
func (double y, void * params)
{
   double temp;
   double * qsat = (double *)params;

   gsl_interp_accel *acc = gsl_interp_accel_alloc ();
   gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, rapmax+1);

   gsl_spline_init (spline, rapgrid, qsat, rapmax+1);

   temp = gsl_spline_eval (spline, y, acc);

   gsl_spline_free (spline);
   gsl_interp_accel_free (acc);

   if ( gsl_isnan(temp) == 1) 
      printf("temp is n.a.n. for y = %f\n",y);
            
   return log(temp*temp);

}



/*********************************************************************
 *    analytic: analytic answer for the  derivative dlog Qs(Y)^2/dY
 *********************************************************************/
double 
analytic (double y, int order)
{
   double gc = 0.6275;
   double xgc = 3.06457;
   double xppgc = 48.5176;
   double pi = 3.14159265359;
   double as = 0.2;
   double sqy=sqrt(y);
   double deriv;
   
   if(order==3) {
      deriv = as * (xgc/gc) - (1.5/gc) / y
                   + 3.0 / (gc*gc) * sqrt(0.5*pi/(as*xppgc)) /(sqy*sqy*sqy);
   }
   else {
      deriv = as * (xgc/gc) - (1.5/gc) / y;
   }
      
   return deriv;
}


/*********************************************************************
 *    reddqs: reduced dlog Qs(Y)^2/dY for FIXED coupling
 *********************************************************************/
double 
reddqs (double y,double dlqs, int order)
{
   double gc = 0.6275;
   double xgc = 3.06457;
   double xppgc = 48.5176;
   double pi = 3.14159265359;
   double as = 0.2;
   double sqy=sqrt(y);
   double result;
   
   if(order==3) {
      result = (dlqs - as * (xgc/gc) + (1.5/gc) / y) * (y*sqy)
             / (1.5 / (gc*gc) * sqrt(2.0*pi/(as*xppgc)));
   }
   else {
      result = 2.0*gc/3.0 * y * (dlqs - as * (xgc/gc));
   }
   return result;
}


/*********************************************************************
 *    reddqs2: reduced dlog Qs(Y)^2/dY for RUNNING coupling
 *********************************************************************/
double 
reddqs2 (double y,double dlqs, int order)
{
   double gc = 0.6275;
   double xgc = 3.06457;
   double xppgc = 48.5176;
   double xi = 2.3381;
   double b = 0.75;
   double sqy=sqrt(y);
   double result;
   
   double c_a = sqrt(2.0*xgc/(b*gc));
   double c_b = 0.75 * xi * exp(log(xppgc/sqrt(2.0*b*gc*xgc))/3.0);
   
   if(order==3) {
      result = 0.0;
   }
   else {
      result = 2.0 * sqy / c_a * (dlqs + c_b/6.0 * exp(-(5.0/6.0)*log(y)));
   }
   return result;
}


/*********************************************************************
 *    bisect:  find the k in [min,max] for which phifunc(k,i) = value
 *             assuming phi is monotonously decreasing
 *********************************************************************/
double
bisect(double value, double min, double max, int i)
{
   double qs;
   double k;
   int comp;
   
   k = exp(0.5*(log(min)+log(max)));
   comp = gsl_fcmp(value,phifunc(k,i),1e-6);
   
   if (comp == 0) {
      qs = k;
   } 
   else if (comp == 1) {
      qs = bisect(value, min, k, i);
   }
   else if (comp == -1) {
      qs = bisect(value, k, max, i);
   }
   else {
      fprintf(stderr,"bisect: something fishy!");
      qs=-1.0;
   }
   return qs;
}


/*********************************************************************
 *    phifunc: return phi(k, i_y) by spline-interpolation in phigrid
 *********************************************************************/
double
phifunc(double k, int iy)
{
   double phi;

   gsl_interp_accel *acc = gsl_interp_accel_alloc ();
   gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, umax+1);

   gsl_spline_init (spline, kgrid, phigrid[iy], umax+1);

   phi = gsl_spline_eval (spline, k, acc);

   gsl_spline_free (spline);
   gsl_interp_accel_free (acc);

   if ( gsl_isnan(phi) == 1) 
      printf("phi is n.a.n. for k = %f\n",k);
            
   return phi;
}   



/*********************************************************************
 *    read_grid: 
 *********************************************************************/
void 
read_grid (char *flname, int *maxy, int *maxu, double *m1, double *m2)
{
   FILE *fp;
   int i=0;
   int y=0;
   int no;
   char dummy[200];
   
   if ( (fp = fopen (flname, "r")) == NULL ) {
	   printf("Error opening file %s.\n", flname);
	   exit(EXIT_FAILURE);
   }
	   
   printf("# Reading data from %s\n", flname);
   (void)fgets(dummy,200,fp);

   while (!feof(fp)) {
      i=0;
      do {
         if((no = fscanf(fp,"%le %le %le %le\n",
         &kgrid[i], &ugrid[i], &rapgrid[y], &phigrid[y][i])) == 4) {
            i++;
         }
         
      } 
      while (!feof(fp) && ugrid[i-1]<=1.0-1.0e-6);

      fscanf (fp, "\n");
      y++;
   }
   fflush(fp);  
   fclose(fp);
   
   *maxy = y-1;
   *maxu =i-1;
   *m1 = -2.0 * log(kgrid[0]);
   *m2 =  2.0 * log(kgrid[i-1]);
   
   return;
   
}
