/************************************************************************* 
 *
 *    bkeqn.c - Main program for numerical solution of the BK equation.
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
#include "parsecfg.h"
#include "bk.h"

void   set_initcond  (double, int);
void   evolve        (void);
void   randomize     (void);

int rapidity = 0;

int    uvMax         = UVMAX;          
int    chMax;        
double lnKm1         = LNKMIN;         
double lnKm2         = LNKMAX;         
int    yMax          = YMAX;           
double deltaY        = YDELTA;         
double odeAccuracy   = ACCURACY;       
double odeStepsize   = STEPSIZE;       
double alphaBar      = ALPHABAR;       
int    alphaRunning  = ALPHARUN;       
int    nFlavour      = NFLAVOUR;       
double lambdaQcd     = LAMBDA;         
double alphaFreeze   = FREEZE;         
int    bkKernel      = BKERNEL;        
int    printSplit    = PRINT_SPLIT;    
int    icChoice      = IC_CHOICE;      
double icPosition    = IC_POSITION;    
char   *icFile;                        
char   *outputFile;                          

int   check_power(int num);
char *xstrdup (const char *str);


/*********************************************************************
 *    main: set up parameters and run
 *********************************************************************/
int main (int argc, char **argv) 
{
   int i;
   char *cfgfilename;
   
   
   /* Define an array of cfgStruct's containing all the 
      variables that can be read from the config-file
      This uses the parsecfg package.  */
   cfgStruct cfg[] = {		
      {"KMAX",          CFG_INT,    &uvMax},                            
	   {"Y_GRID_MAX",    CFG_INT,    &yMax},                            
	   {"N_FLAVOUR",     CFG_INT,    &nFlavour},                            
	   {"BK_KERNEL",     CFG_INT,    &bkKernel},                            
	   {"IC_CHOICE",     CFG_INT,    &icChoice},                           
	   {"ALPHA_RUNNING", CFG_BOOL,   &alphaRunning},                     
	   {"PRINT_SPLIT",   CFG_BOOL,   &printSplit},                            
	   {"LN_K_MIN",      CFG_DOUBLE, &lnKm1},                            
	   {"LN_K_MAX",      CFG_DOUBLE, &lnKm2},                            
	   {"DELTA_Y",       CFG_DOUBLE, &deltaY},                            
	   {"ODE_ACCURACY",  CFG_DOUBLE, &odeAccuracy},                            
	   {"ODE_STEPSIZE",  CFG_DOUBLE, &odeStepsize},                            
	   {"ALPHA_BAR",     CFG_DOUBLE, &alphaBar},                            
	   {"LAMBDA_QCD",    CFG_DOUBLE, &lambdaQcd},                            
	   {"ALPHA_FREEZE",  CFG_DOUBLE, &alphaFreeze},                            
	   {"IC_POSITION",   CFG_DOUBLE, &icPosition},                            
      {"IC_FILE",       CFG_STRING, &icFile},                        
      {"OUTPUT_FILE",   CFG_STRING, &outputFile},                        
	   {NULL,            CFG_END,    NULL}	   
   };


   /* If no config filename given on command line, use default settings */
   if( argc == 1 ) {
      printf("-----------------------------\n");            
      printf("Note: using default settings!\n");            
      printf("-----------------------------\n");            
      uvMax = UVMAX;       
      yMax = YMAX;        
      nFlavour = NFLAVOUR;    
      bkKernel = BKERNEL;    
      icChoice = IC_CHOICE;   
      alphaRunning = ALPHARUN;
      printSplit = PRINT_SPLIT;  
      lnKm1 = LNKMIN;       
      lnKm2 = LNKMAX;       
      deltaY = YDELTA;      
      odeAccuracy = ACCURACY; 
      odeStepsize = STEPSIZE; 
      alphaBar = ALPHABAR;   
      lambdaQcd = LAMBDA;   
      alphaFreeze = FREEZE;  
      icPosition = IC_POSITION;         
   }
   
   
   /* Else if filename given, read settings from the config-file  */
   else if( argc == 2 ) {
      if ((cfgfilename = xstrdup(argv[1])) == NULL ) {
         fprintf(stderr,"Error: xstrdup failed in main()!\n");
         exit(EXIT_FAILURE);
      }
       
      for(i=1 ; i<(int)strlen(cfgfilename)+26 ; i++) printf("-");
      printf("\n");            
      printf("Note: using config file %s!\n",cfgfilename);            
      for(i=1 ; i<(int)strlen(cfgfilename)+26 ; i++) printf("-");
      printf("\n");            
     
      if (cfgParse(cfgfilename, cfg, CFG_SIMPLE) == -1) {
         fprintf(stderr,"Error reading configuration file!\n");
         exit(EXIT_FAILURE);
      }

      free(cfgfilename);
   }
   
   
   /* Else too many args */
   else {
      fprintf(stderr,"Usage: %s <optional config-filename>\n",argv[0]);
      exit(EXIT_FAILURE);
   }   
   
     
      
   /* Set filenames if not already set */
   if( icFile == NULL ) {
      icFile = xstrdup(IC_FILE); 
   } 
   if( outputFile == NULL ) {
      outputFile = xstrdup(OUTPUT_FILE); 
   } 
   if( icFile == NULL || outputFile == NULL ) {
      fprintf(stderr,"Error: xstrdup failed in main()!\n");
      exit(EXIT_FAILURE);
   } 


   if( check_power(uvMax) != 0 ) {
      fprintf(stderr,"Error: KMAX must be a power of two!\n");
      exit(EXIT_FAILURE);
   }

   chMax = 2*uvMax;

 
   printf("uvMax = %d\n",uvMax);            
   printf("chMax = %d\n",chMax);        
   printf("lnKm1 = %g\n",lnKm1);            
   printf("lnKm2 = %g\n",lnKm2);            
   printf("yMax = %d\n",yMax);        
   printf("deltaY = %g\n",deltaY);           
   printf("odeAccuracy = %g\n",odeAccuracy);   
   printf("odeStepsize = %g\n",odeStepsize);   
   printf("alphaBar = %g\n",alphaBar);        
   printf("alphaRunning = %d\n",alphaRunning);  
   printf("nFlavour = %d\n",nFlavour);        
   printf("lambdaQcd = %g\n",lambdaQcd);      
   printf("alphaFreeze = %g\n",alphaFreeze);   
   printf("bkKernel = %d\n",bkKernel);        
   printf("printSplit = %d\n",printSplit);    
   printf("icChoice = %d\n",icChoice);        
   printf("icPosition = %g\n",icPosition);    
   printf("icFile = %s\n",icFile);                      
   printf("outputFile = %s\n",outputFile);                
   
   
   /* Boring setup complete, now run the solver */
   set_initcond(icPosition, icChoice);
   evolve();
   
   return 0;
}
