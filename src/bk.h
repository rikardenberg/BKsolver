/************************************************************************* 
 *
 *    bk.h - Define default settings and declare external variables.
 *           Don't use this to change values, use config-file instead.
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

/* Basename for output file */
#define OUTPUT_FILE "bkresult"

/* Number of u-intervals. MUST be a power of 2 !  */
#define UVMAX 256   

/* Max and -min of log(k^2) */
#define LNKMAX  138.0
#define LNKMIN  20.0

/* Number of Y-intervals */
#define YMAX 10

/* Length of each Y-interval (note that this is not the ODE step size!) */
#define YDELTA 5

/* Accuracy and initial step-size in ODE solution */
#define ACCURACY 1e-3
#define STEPSIZE 1e-3

/* Alpha-bar value of fixed coupling */
#define ALPHABAR 0.2

/* Alpha-bar running (1) or fixed (0) ? */
#define ALPHARUN 0

/* Running coupling-related */
#define NFLAVOUR 3
#define LAMBDA 0.2
#define FREEZE 0.5

/* Which kernel - BFKL: 0 or Two-pole: 1 or BFKL with approx CC: 2*/
#define BKERNEL 0

/* Choice of initial condition
    Read from file:   0
    Step function:    1
    Smooth step:      2
    Gaussian:         3 
    1/sqrt(k):        4    (which is not-steep-enough...)  */
#define IC_CHOICE 0

/* File containing table of initial condition if IC_CHOICE = 0 */
#define IC_FILE "mv-k.front.dat"

/* Position of initial condition if IC_CHOICE > 0 */
#define IC_POSITION 1.0

/* Print data to single file (0) or multiple files (1) */
#define PRINT_SPLIT 0

/* External variables */
extern double  **phigrid;
extern double  *rapgrid;
extern double  *ugrid;
extern double  *ws;
extern double  *abar;
extern int     rapidity;

extern int    uvMax;          
extern int    chMax;       
extern double lnKm1;          
extern double lnKm2;          
extern int    yMax;       
extern double deltaY;         
extern double odeAccuracy;    
extern double odeStepsize;    
extern double alphaBar;       
extern int    alphaRunning;   
extern int    nFlavour;       
extern double lambdaQcd;      
extern double alphaFreeze;    
extern int    bkKernel;       
extern int    printSplit;     
extern int    icChoice;       
extern double icPosition;     
extern char   *icFile;       
extern char   *outputFile;                    

