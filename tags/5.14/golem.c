/************************************************************************
 *                                                                            
 *  GOLEM
 *
 *  Version 5.14
 *  September, 1999
 *
 *  NOTE: This code is freely available for academic research and
 *        teaching purposes only. Commercial use of this program in
 *        whole or in part is expressly prohibited without written
 *        consent from the author. For further information, contact
 *        Greg Tucker, Dept. of Civil & Environmental Engineering,
 *        Mass. Inst. of Technology, Cambridge, MA 02139, email
 *        gtucker@mit.edu.
 *
 *  Copyright (C) 1993-2000 Gregory E. Tucker
 *
 *  Additions and changes from 5.13 to 5.14:
 *    Incorporated the following updates accidentally made to older
 *       version (5.12):
 *    1 - 3/98: fixed an obscure bug in slopeCollapse: if collapse occurs
 *             any ponding info is lost when drainage dirs are reset.
 *             This led to problems when both br and fluvial were "on".
 *             A further problem was that downstream smoothing was
 *             being applied even to flooded cells, with the result that
 *             a positive sediment flux was imparted to flooded cells.
 *    2 - 3/98: modified DiffuseExplIrreg so that it can operate in
 *             either "weathering-limited" or "diffuseall" mode, instead
 *             of only the former.
 *    Other fixes:
 *    3 - 6/99: fixed bug in DiffuseAll: dhin array needed to be
 *             initialized over ALL entries, not just 1:XGRID (etc)
 *    4 - 9/99: fixed bug in DiffuseExplIrreg caused by 2 (above):
 *             negative sed thickness caused by large divergence and
 *             small soil thickness.
 *    5 - 9/99: removed comment accidentally left behind after debugging
 *             #4!
 *
 *  Changes to 5.14:
 *    1 - 6/00: d'oh! accidentally deleted call to SlopeCollapse() during
 *             the last debugging session...
 *    2 - 10/00: minor change to ThresholdTransport to improve efficiency.
 *    3 - 3/01: added a few comments and removed an unused vbl in
 *              checkSlopes() function
 *    4 - 3/02: fixed bug in EstimateDtMax: dimension of dht[][] incorrect.
 *              Changed from dht[XGRID+2][YGRIDP1] to dht[XGRID+2][YGRID+2]
 *
 *  $Id: golem.c,v 1.5 2002/03/29 11:30:59 gtucker Exp $
 *
 ************************************************************************/


                /* HEADER FILES */

#include <stdio.h>  /* For file I/O */
#include <stdlib.h> /* For qsort routine and random # fn's */
#include <malloc.h> /* For allocating lithology array and pond lists */
#include <math.h>   /* Power and trig functions, etc */
#include <string.h>


                /* DEFINED CONSTANTS */

#define VERSION "Version 5.14 with changes 1-4, Mar 2002"
#define XGRID 62        /* Grid size in X, excluding boundaries */
#define YGRID 62        /* Grid size in Y, excluding boundaries */
#define XGRIDP1 63       /* XGRID + 1 (= coord of east boundary) */
#define YGRIDP1 63      /* YGRID + 1 (= coord of north boundary) */
#define GRIDSIZE 3844    /* Number of cells, excluding boundaries */
#define TRUE 1          
#define FALSE 0
#define NEWPOND 1       /* Code used in ponding algo */
#define OLDPOND 2       /* Code used in ponding algo */
#define BASIN 3         /* Code for cell identified as closed basin */
#define PLATEAU 0       /* Code for uniform uplift, two fixed boundaries */
#define BLOCK 1         /* Code for uniform uplift, one fixed boundary */
#define TILTBLOCK 2     /* Code for tilt block uplift, one fixed boundary */
#define DYNAMIC 3       /* Code for dynamic uplift, ie uplift = avg denud */
#define ERRFUNC 4       /* Code for error function uplift */
#define HALFDOME 5      /* Code for half dome uplift */
#define ERFSUB 6        /* Code for error function time-dep subsidence */ 
#define FOLDS 7         /* Code for sinusoidal folding */
#define STRIKESLIP 8    /* Code for strike-slip motion */
#define SECPERHOUR 3600.0       /* No. of seconds in an hour */
#define HOURPERYEAR 8760.0      /* No. of hours in a year */
#define CHANNEL 1       /* Code for a channel cell */
#define HILLSLOPE 0     /* Code for a hillslope cell */
#define BELOWTHRESH 1   /* Code for pixel below erosion threshold */
#define ALLUVIAL 2      /* Code for an alluvial channel cell */
#define CAPLIM 3        /* Code for a capacity-limited bedrock channel cell */
#define BEDROCK 4       /* Code for a bedrock channel cell */
#define RHO 1000.0      /* Density of water (kg/m^3) */
#define FIXED_ZERO 0    /* East bdy code: fixed at zero elevation */
#define FIXED_THICKNESS 1   /* East bdy code: fixed rock/sed thickness */
#define NO_FLUX 2       /* East bdy code: no flux (zero slope maintained) */ 
#define SLOPE_EXTRAP 3  /* East bdy code: (see adjustEastBounday() ) */
#define GRAV 9.81       /* Gravitational acceleration, m/s^2 */
#define SOIL 1          /* Code for a soil-landslide (see slopeCollapse() ) */
#define ROCK 2          /* Code for a rock-landslide (see slopeCollapse() ) */
#define VERYSMALL 0.0001
#define INIT 1          /* Code for initialization, used by foldUplift() */
#define ROOT2 1.4142136 /* The square root of two */
#define NumUpliftTypes 8    /* Number of valid uplift types */
#define NumEastBdyOpts 4    /* Number of options for east boundary cond'n. */
#define AverageElev 0       /* Flag value: track avg cell elevation */
#define ChannelElev 1       /* Flag value: track channel elevations */
#define CatchmentMode      0   /* Model runs in catchment-scale mode */
#define RegionalScaleMode  1   /* Model runs in regional scale mode */
#define NoFluxBoundary  1   /* Code for no-flux boundary cell, if irreg. bdy*/
#define OpenBoundary    2   /* Code for a legal exit point */
#define BAGNOLD         1   /* Code for Bagnold bedload transp. function */
#define STREAMPWR       2   /* Code for simple strm. pwr. transp. function */
#define THRESHOLD       3   /* Code for generic threshold power transport fn*/
#define MPM             4   /* Code for Meyer-Peter Mueller transport fn */
#define WEATHLIMDIFF    1   /* Code for weathering-limited diffusion */
#define STEPWAVE 2          /* Code for step-function variation in precip */
#define PI 3.1415926        /* Approx value of pi */ 
#define TWOPI 6.2832        /* Approx value of 2*pi */ 


/*                          GLOBAL CELL VARIABLES 
 * These arrays store variables that are tracked for each grid cell location,
 * such as elevation, flow direction, sediment cover thickness, etc.
 * For the variables that are declared as array pointers, such as
 * **tau, memory is dynamically allocated only when needed.
 */ 
double elev[XGRID+2][YGRID+2];          /* cell elevations */
double load[XGRID+2][YGRID+2];          /* sediment infux from upstream */
double slope[XGRID+2][YGRID+2];         /* steepest slope */
double disch[XGRID+2][YGRID+2];         /* water discharge */
double delx[XGRID+2][YGRID+2];          /* stream channel length betwn cells*/
int nbrx[XGRID+2][YGRID+2];             /* x-coord of lowest neighbor cell */
int nbry[XGRID+2][YGRID+2];             /* y-coord of lowest neighbor cell */
int area[XGRID+2][YGRID+2];             /* Drainage area in cells */
double chansed[XGRID+2][YGRID+2];       /* Sediment cover thickness */
float sla[XGRID+2][YGRID+2];            /* Sed thickness: 3 pt moving avg */
float slb[XGRID+2][YGRID+2];            /* Sed thickness: 3 pt moving avg */
int rockType[XGRID+2][YGRID+2];         /* Lithology of uppermost rock layer*/
double rock[XGRID+2][YGRID+2];          /* Thickness of uppermost rock layer*/
double base[XGRID+2][YGRID+2];          /* Depth (neg) to base of rock col. */
char flood[XGRID+2][YGRID+2];           /* Flooding status indicator */
float delz[XGRID+2][YGRID+2][3];        /* 3-pt moving average denud rate */
char chan[XGRID+2][YGRID+2];            /* Channel (or hillslope) type */
char boundary[XGRID+2][YGRID+2];        /* Map of active portion of grid */
double chanDepth[XGRID+2][YGRID+2];     /* Channel depth */
char **slp_fail;                        /* Slope failure locations */
double **tau;                   /* Bed shear stress: used only when the
                                   Bagnold transport function is used */
double **hillslopeReservoir;    /* Volume of rock and sediment
                                   in each cell apart from the channel
                                   sediment and rock reservoirs.
                                   Expressed in units of cell elevation
                                   (i.e. volume required to fill a
                                   cell to X meters depth) */
double **meanElev;              /* Mean elevation of cells. Used in
                                   regional-scale mode when computing
                                   channel bed elevations rather than
                                   mean elevations. */


        /* STRUCTURES AND TYPE DEFINITIONS */

struct CellCoord        /* This structure stores the coords of a cell */
{
    int x;
    int y;
};

struct PondCell         /* Element of linked list used in ponding algo */
{
    int x;
    int y;
    struct PondCell *next;
};

struct CellPoint        /* Coordinates AND elevation of a cell */
{
    int x;
    int y;
    double z;
};

/*
 *  Lithology structure: Version 5.10 adds fields, Kmn and tc. Kmn is the mean,
 *    unmodified K value---that is, unmodified by changes in precip frequency.
 *    If precip frequency/duration varies, the "effective" K value will also
 *    vary. Thus, Kmn stores the mean, and K the value at any particular
 *    time, as computed by the CalcPrecip function. If precip doesn't vary
 *    in time, K and Kmn will be identical. tc is a bedrock erosion threshold,
 *    the units of which depend on exponents mb and nb. It represents that
 *    value of Q^mb S^nb at which erosion just begins to be initiated.
 */  
struct Lithology        /* Lithology record */
{
    double Kmn;      /* Erodability by flowing water (mean) */
    double K;        /* Erodability by flowing water (time-varying) */
    double T;        /* Initial thickness */
    double W;        /* Bedrock weathering rate */
    double Scr;      /* Failure slope */
    double tc;       /* Erosion threshold */
} *layer; 

typedef char *String;


        /* OTHER GLOBAL VARIABLES */

int Dt,             /* Time step size (years) */
    writeInterval,  /* Output interval in time-steps */
    start_time=0,   /* Starting time index (>0 for restart runs) */
    end_time,       /* Ending time index (=total #ts plus start time) */
    uptime,         /* Duration of tectonic uplift   */
    uptype,         /* Type of tectonic uplift       */
    east_boundary,  /* East (lefthand) boundary condition */
    n_layers,       /* Number of rock layers */
    numActiveCells, /* Number of active (non-boundary) grid cells */
    streamCaptureOccurred=TRUE; /* Flag indicating whether a change in 
                       drainage direction has occurred during the current
                       time step. */

int pond,           /* Flag for ponding option: 0=no ponding */
    optIrregularBoundary,   /* Flag indicating irregular boundary */
    optLayeredLithology=1,  /* Flag indicating layered lithologies */
    optTransFn,     /* Code for type of fluvial transport function */
    optVarPrecip,   /* Option for time-varying precip intensity */
    nPtSources,     /* No point sources for discharge @ east bnd */
    *inputLoc,      /* Location(s) of discharge point sources @ east bnd */
    *Ainput,        /* Drainage area feeding point source inputs */
    permChannelA,   /* Minimum drainage area for permanent, non-weathering
                       channel (in cells). Weathering is not applied to cells
                       with area above this value. */
    fault,          /* x-coordinate of fault location, used by uplift fn's */
    erfoffset,      /* Offset for "error function" uplift, in cells */
    tsubstart,      /* Start time of subsidence, used in erfSubsidence */
    optMinSlope,    /* This option keeps elev's above 0 in bedrk-only mode */
    optDiffuse,     /* Option for diffusion function type */
    optAreaDepLS,   /* Option for area-dependent landsliding */
    optSat=0,       /* Option for saturation-excess runoff */
    ngrabens,       /* Number of half grabens at boundary */
    grabenLength;   /* Length in cells of half grabens near boundary */

double  kf,         /* Fluvial sediment transport coefficient */
    kfmean,         /* Mean transport coef (varies w/ precip intensity to
                       conserve total precip---reflects storm recurrence) */
    gkt,            /* Shear stress coefficient used in ThresholdTransport */
    mf, nf,         /* Discharge and slope exp's for fluvial sed transport */
    gpf,            /* Exponent for excess shear stress, ThresholdTransport */
    gtaucrit,       /* Crit shear stress for sediment, ThresholdTransport */
    gkc,            /* Chan width coef: MPMTransport */
    gkdp,           /* Chan depth coef */
    gmdp,           /* Chan depth exponent */
    kbeff,          /* "Effective" bedrock coef (varies w/ precip frequency)*/
    mb, nb,         /* Discharge and slope exp's for bedrock erosion */
    kd,             /* Hillslope transport (diffusivity) coefficient */
    upliftRate,     /* Uplift rate constant (m/yr) */
    denud,          /* Used to calculate average denudation rate */
    fault_to_pivot, /* Distance from fault to tiltblock pivot, m */
    vf,             /* Velocity of westward fault zone migration */
    erfdist,        /* Distance of uplift zone for "error function" uplift */
    chanThreshold,  /* Threshold of channel initiation */
    mci,            /* Channel initiation fn: Q^mci */
    nci,            /* Channel initiation fn: S^nci */
    mw,             /* Depth at which weathering rate decays to 1/e */
    soil_slope,     /* Failure slope for soil layer (for area-dependent
                       landsliding, this is the angle of internal friction */
    angleOfRepose,  /* Angle of repose for material released in a landslide */
    kls,            /* Coefficient for computing critical area in area-dep't
                       landsliding */
    precip_flow,    /* Volumetric flow rate from precip on 1 cell (m3/yr) */
    precip_mean,    /* Mean precip on 1 cell (used if time-varying) */
    precip_var,     /* Maximum variation of precip intensity */
    gp_mean_ann,    /* Mean annual flow/area, used by MPMTransport */
    lambdap,        /* Time scale of precip variation (one period, in Dt's) */
    minS,           /* Minimum req'd channel gradient for br chan (optional) */
    dt_min,         /* Minimum time slice (yrs) into which a time step 
                       can be broken */
    dtMinApplied,   /* Smallest sub-time slice actually used during a 
                       given time step (used only in tracking) */
    downstrmWt,     /* Weighting factor for solution averaging in 
                       EstimateDtMax (improves fluvial transport stability) */
    upstrmWt,       /* Weighting factor for solution averaging in 
                       EstimateDtMax (improves fluvial transport stability) */
    soilTrans,      /* soil transmissivity (m2/yr) (for saturation runoff) */
    gwrchg,         /* groundwater recharge rate (for sat runoff) */
    sliprate,       /* slip rate in m/yr for StrikeSlip "uplift" function */
    dx,             /* Horizontal grid spacing, in meters */
    dxx;            /* Diagonal distance between cells */

double (*Transport)();  /* pointer to sed transport function of choice */

struct CellCoord cell[GRIDSIZE];  /* List of cell coordinates */

double cellarea,    /* Area of a cell = dx*dx */
    chanarea,       /* Area of a cell subject to active channel erosion and
                       deposition, m^2 */
    R,              /* Coefficient for ADI diffusion solver (UNUSED) */
    *uca,*ucb,*ucc, /* Tridiagonal matrix coef's for ADI solver  (UNUSED)*/
    *vca,*vcb,*vcc; /*     "          "        "  "   "    "     (UNUSED)*/

int elevToTrack,    /* Flag indicating whether ELEV and ROCK refer to
                       the average cell elevation (AverageElev) or to
                       the principal channel elevation (ChannelElev) */
    slopeCollapseActive,    /* Flag indicates whether to bother with
                       computing slope collapse */
    mode,           /* Simulation mode: CatchmentMode (channels and
                       hillslopes resolved) or RegionalScaleMode (channels
                       modeled and hillslopes parameterized) */
    numChannelCells=GRIDSIZE;   /* Number of cells that contain channels */


        /* VARIABLES USED IN 2-D ELASTIC BEAM FLEXURE CALCULATION */

int flexure;        /* Flag indicating activation of flexure */
double flexr,       /* Flexural rigidity */
    alpha,          /* Flexural parameter */
    rhom,           /* Mantle density */
    rhoc,           /* Crustal density (uniform for all lithologies) */
    rhos,           /* Regolith density */
    w[XGRID+2],     /* Deflection (m), positive **? */
    V[XGRID+2],     /* Crustal load (N/m) */
    plat_w[XGRID+2],/* Deflection due to uneroded plateau beyond east bdy */
    plat_init[XGRID+2],   /* Initial deflection from plateau */
    shelf_w[XGRID+2],     /* Deflection due to shelf sediment loading */
    w_init[XGRID+2],      /* Initial deflection */
    init_up[XGRID+2],     /* Initial uplift */
    plat_thick,     /* Thickness of interior plateau */
    unitload,       /* = alpha^3/8D */
    sed_yield,      /* Sediment outflux at west boundary. */
    shelf_param,    /* Scaling factor converts sed yield to shelf load */
    shelf_load,     /* Sediment load on shelf (in each 5km cell) */
    sselev;         /* Elevation of unloaded, unflexed plate above SL */


        /* OUTPUT FILE POINTERS */

FILE *elevfp,       /* Cell elevations for each recorded time step */
    *mnelevfp,      /* Mean elevation (used in RegionalScaleMode) */
    *netfp,         /* Drainage net info: channel flags & neighbr coords */
    *solfp,         /* Soil thickness */
    *afp,           /* Drainage area */
    *qfp,           /* Discharge (only written if optSat) */
    *sedfp,         /* Sediment yield at west boundary */
    *denudfp,       /* Average denudation rate */
    *lithfp,        /* Surface lithologies */
    *chanfp,        /* Channel types (bedrock vs alluvial) */
    *timefp,        /* Current time step */
    *flexfp,        /* Flexural profile */
    *slpfp,         /* Regions of slope failure */
    *taufp,         /* Bed shear stresses */
    *lyrfp;         /* Thickness of uppermost rock layer */
char restart_name[50];
char timefile[50];


        /* FUNCTION DECLARATIONS */

void ReadLithology();
void FlexInit();
void ResetVariables();
void ShiftAverages();
void adjustEastBoundary();
void FindNextCell();
double GenericTransport();
double ThresholdTransport();
double MPMTransport();
double Bagnold();
double ErodeBedrock();
void DiffuseADI();            /* ADI solver; not currently in use */
void DiffuseExpl();
void DiffuseExplIrreg();
void DiffuseAll();
void wrapNSBoundaries();
struct PondCell *AddList();
struct Lithology *LithVector();
int CompareAreas();
int CompareChannelAreas();
static void tridag();
void Flex();
void SlopeCollapse();
void foldUplift();
void problem();
char *date();
void WriteRestartFile();
void ADIinit();      /* Initialization for ADI solver; not currently in use */
double *dvector();
double **dmatrix();
char **cmatrix();
void StrikeSlip();
void shunt(int k,int m,struct CellCoord temp,struct CellCoord array[]);
void heapsort(int n,struct CellCoord array[]);


/***********************************************************************
 *      MAIN FUNCTION                                                   
 ***********************************************************************
 *                                                                      
 *  main:       The main function reads in initial data, calls the      
 *              initialize function, and implements the main            
 *              time loop.
 *
 *      Arguments: argc & argv -- standard C command-line argument
 *              variables.
 *                                                                      
 *      Calls: problem, initialize, Flex, FindNextCell, DoPonds,
 *              FindChannels, outElevData, outMeanElevData,
 *              outNetData, outSedData, outChanData, outDelzData,
 *              outSlpData, outLithologyData, outLayer1Data,
 *              outQData, outFlexData, qsort [library function],
 *              Weather, DiffuseExplIrreg, DiffuseExpl, 
 *              SlopeCollapse, Fluvial, Uplift, WriteRestartFile
 *                                                                      
 ************************************************************************/

int main(argc,argv)
unsigned int argc;      /* Number of arguments given on command line */
char **argv;            /* Command line arguments (array of strings) */
{
    int i, j;           /* Counters for x and y coordinates */
    int time,           /* Current time step */ 
        silent=FALSE,   /* Flag for tty output */
        next_write;     /* Next time step for output to data files */

    /* SAY HELLO */
    printf( "\n******** G O L E M ********\n\n" );
    printf( "Written by Greg Tucker\n" );
    printf( "Copyright (C) 1996-2002, G.E. Tucker\n\n" );
    printf( "%s\n", VERSION );
    printf( "Compiled for %d by %d grid (including boundaries)\n\n", 
	    XGRID+2, YGRID+2 );

    /* MAKE SURE AN INPUT FILE HAS BEEN SPECIFIED */
    if( argc < 2 ) problem( "No input file specified." );

    /* CHECK FOR "SILENT" OPTION */
    if( argc > 2 && argv[2][1]=='s' ) silent=TRUE;

    /* INITIALIZATION */
    initialize( argv[1] );
    next_write = (start_time==0) ? 0 : start_time + writeInterval; 

    /* MAIN LOOP */
    for( time=start_time; time<=end_time; time++ ) 
    {
        /* Report the current time step to screen and/or file */
        if( !silent ) fprintf( stderr, " %d\n", time ); 
        if( fmod((float)time, 100.0)<1 ) {
            timefp=fopen(timefile,"wt");
            fprintf( timefp, "Time = %d\n",time); 
            fclose(timefp);
        }

        /* Calculate flexural isostasy, if applicable */
        if( flexure ) Flex();

        /* Reset variables for new time-step and compute new elevations */
        ResetVariables();

        /* Compute flow directions and if necessary drain closed basins */
        for( i=1; i<=XGRID; i++ )
            for( j=1; j<=YGRID; j++ )
                if( !boundary[i][j] ) FindNextCell( i,j );
        if( pond && streamCaptureOccurred )
            for( i=1; i<=XGRID; i++ )
                for( j=1; j<=YGRID; j++ )
                    if( flood[i][j]==BASIN && !boundary[i][j] )
                        DoPonds( i, j );

        /* If any flow dir's have changed, compute new drainage patterns */
        /* Also locate channels if in CatchmentMode */
        if( streamCaptureOccurred ) StreamTrace();
        CalcDischarge();
        if( mode==CatchmentMode && chanThreshold>0 ) FindChannels();

        /* Write output to files at selected intervals */
        if( time==next_write )
        {
            outElevData( time );
            if( mode==RegionalScaleMode && elevToTrack==ChannelElev )
                outMeanElevData( time );
            outNetData( time );
            if( time > 0 ) 
            {
                outSoilData( time );
                outSedData( time );
                outChanData( time );
                outDelzData( time );
                if( slopeCollapseActive ) outSlpData( time );
            }
            if( n_layers > 1 ) outLithologyData( time );
            if( n_layers > 1 ) outLayer1Data( time );
            outQData( time );
            if( flexure ) outFlexData( time );
            if( optTransFn==BAGNOLD || optTransFn==THRESHOLD 
                || optTransFn==MPM) 
                outTauData( time );
            next_write += writeInterval;
        }
        ShiftAverages();        /* Used for moving avg sed & denud rate */

        /* Sort the cells by area and compute hillslope processes */
        if( mode==CatchmentMode )
        {
	  if( streamCaptureOccurred ) /* Only re-sort when necessary */
	    /*printf("List before sorting:\n");
	    for( i=0; i<numActiveCells; i++ )
	      printf("%d %d %d\n",cell[i].x,cell[i].y,area[cell[i].x][cell[i].y]);
	      heapsort( GRIDSIZE, cell );
	    printf("List after sorting:\n");
	    for( i=0; i<numActiveCells; i++ )
	    printf("%d %d %d\n",cell[i].x,cell[i].y,area[cell[i].x][cell[i].y]);*/
	    qsort( cell, numActiveCells, sizeof( cell[0] ), 
	      CompareChannelAreas );
            Weather();
            if( kd > 0.0 ) {
                if( optIrregularBoundary ) DiffuseExplIrreg();
                else if( optDiffuse==WEATHLIMDIFF )
                    DiffuseExpl();
                else DiffuseAll();
            }
            if( slopeCollapseActive ) SlopeCollapse();
        }
        else
        {
            if( streamCaptureOccurred ) /* Only re-sort when necessary */
            {
	      heapsort( numActiveCells, cell );
	      /*if( optIrregularBoundary )
                    qsort( cell, numActiveCells, sizeof( cell[0] ), 
                           CompareChannelAreas );
                else
                    qsort( cell, numActiveCells, sizeof( cell[0] ), 
		    CompareAreas );*/
            }
            /*SedFluxToChannels(); NOT YET IMPLEMENTED! */
        }

        /* Compute time-varying precip intensity, if applicable */
        /*if( optVarPrecip ) CalcPrecip( time );*/

        /* Compute fluvial transport and erosion/deposition */
        Fluvial();

        /* Tectonic uplift, anyone? */
        if( time < uptime ) Uplift( uptype, time );

        streamCaptureOccurred = FALSE;  /* Re-set for the new time-step */

    } /* End of time-step loop */

    WriteRestartFile( time-1 );

    return 0;

} /**** end of main *****************************************************/



/***********************************************************************
 *
 *  ReadInt:    ReadInt() reads a line of text followed by an integer
 *              from the file pointed to by fp. It returns the value
 *              of the integer (the string is discarded).
 *
 *      Argument: fp -- pointer to input file being read.
 *      Returns: the integer read from the file.
 *      Calls: (library functions only)
 *      Called by:  initialize(), foldUplift()
 *                                                                      
 ************************************************************************/

int ReadInt( fp, itemCode, nchars )
FILE *fp;
char *itemCode;
int nchars;
{
    int c, value;
    char headerLine[80];

    fgets( headerLine, 80, fp );
    fscanf( fp, "%d", &value );
    fgetc( fp );
    for( c=0; c<nchars; c++ )
        if( headerLine[c] != itemCode[c] )
        {
            printf("I expected to read the parameter %s, but instead found:\n",
                     itemCode );
            printf("    %s    %d\n",headerLine,value );
            problem("Input file format incorrect\n"); 
        }
    return( value );
}



/***********************************************************************
 *
 *  ReadDouble: ReadDouble() reads a line of text followed by a double
 *              from the file pointed to by fp. It returns the value
 *              of the double (the string is discarded).
 *
 *      Arguments: fp -- pointer to the input file being read.
 *      Returns: the value read from the file.
 *      Calls: (library functions only)
 *      Called by:  initialize()
 *                                                                      
 ************************************************************************/

double ReadDouble( fp, itemCode, nchars )
FILE *fp;
char *itemCode;
int nchars;
{
    float value;
    char headerLine[80];
    int c;

    fgets( headerLine, 80, fp );
    fscanf( fp, "%f", &value );
    fgetc( fp );
    for( c=0; c<nchars; c++ )
        if( headerLine[c] != itemCode[c] )
        {
            printf("I expected to read the parameter %s, but instead found:\n",
                     itemCode );
            printf("    %s    %f\n",headerLine,value );
            problem("Input file format incorrect\n");
        }
    return( (double)value );
}



/***********************************************************************
 *
 *  ReadFloat:  ReadFloat() reads a line of text followed by a float
 *              from the file pointed to by fp. It returns the value
 *              of the float (the string is discarded).
 *
 *      Arguments: fp -- pointer to the input file being read.
 *      Returns: the value read from the file.
 *      Calls: (library functions only)
 *      Called by:  foldUplift()
 *                                                                      
 ************************************************************************/

float ReadFloat( fp, itemCode, nchars )
FILE *fp;
char *itemCode;
int nchars;
{
    float value;
    int c;
    char headerLine[80];

    fgets( headerLine, 80, fp );
    fscanf( fp, "%f", &value );
    fgetc( fp );
    for( c=0; c<nchars; c++ )
        if( headerLine[c] != itemCode[c] )
        {
            printf("I expected to read the parameter %s, but instead found:\n",
                     itemCode );
            printf("    %s    %f\n",headerLine,value );
            problem("Input file format incorrect\n");
        }
    return( value );
}



/***********************************************************************
 *
 *  CompareStrings: This function compares the first three characters
 *                  of two strings. If they're identical, it returns
 *                  a value of TRUE, otherwise FALSE.
 *              
 *      Arguments: the two strings to compare.
 *      Returns: TRUE if the first 3 characters are the same, FALSE
 *               otherwise.
 *      Calls: (none)
 *      Called by:  initialize()
 *                                                                      
 ************************************************************************/

int CompareStrings( str1, str2 )
String str1, str2;
{
    int i, str1start=0;

    /* Get to the first non-blank char in str1 */
    while( str1[str1start]==' ' || str1[str1start]=='\t' ) str1start++;

    /* Compare the 1st 3 characters: if one doesn't match, 
       return FALSE, otherwise TRUE. */
    for( i=0; i<=2; i++ )
        if( str1[i+str1start]!=str2[i] ) return( FALSE );
    return( TRUE );
}



/***********************************************************************
 *
 *  initialize: The initialization does the sort of things you'd expect
 *              it to: it reads in the main input file, including the
 *              names of the stratigraphy and base files, and if needed
 *              the boundary file name as well. The boundary type and
 *              uplift type are read as strings, so they're converted
 *              to number codes using the upliftName[] and
 *              eastBndCond[] String arrays. After reading in stuff from
 *              the main input file, this function then:
 *
 *               o Writes the date, time & code version onto the end of
 *                  the input file
 *               o allocates space for additional arrays if needed
 *               o calls ReadLithology() to set up initial topography
 *                  and stratigraphy
 *               o if needed, reads a boundary map file
 *               o sets up the boundary codes
 *               o initializes a number of variables
 *               o calls functions to set up a list of grid cells,
 *                  open and initialize the output files, and if
 *                  needed initialize the flexure equations
 *
 *      Arguments: infile is the name of the input file.
 *      
 *      Calls:  problem(), ReadInt(), ReadDouble(), CompareStrings(),
 *              foldUplift(), cmatrix(), dmatrix(), ReadLithology(),
 *              MakeSortList(), OpenOutputFiles(), FlexInit()
 *
 *      Called by:  main()
 *                                                                      
 ************************************************************************/

initialize( infile )
char *infile;           /* Name of input file */
{
    int k, i, j,        /* Loop counters */
        restart;        /* Flag indicating new run from restart file */
    int seedValue;      /* Seed for rand # generator, used in foldUplift */
    double basin_width; /* Width of offshore basin if flexure is on */
    double kci;         /* Coef. for chan init fn (eg, constant relating
                           Q and S to shear stress, etc) */
    FILE *ifp;          /* -> the input file */
    float solinit;      /* Initial soil thickness */
    char tmpCharArray[YGRID+2][XGRID+2];  /* Used for binary read of 
                                             boundary data (if needed) */
    char stratfile[50],     /* Stratigraphy/lithology input file name */
        basefile[50],       /* Basal surface input file name */
        outputfile[50],     /* Base name for output files */
        tempstr[80],        /* Temporary string used in reading... note that
                               to read a string after a #, need to read 2
                               strings and discard the 1st */
        boundaryfile[80],   /* File containing irregular boundary map */
        restartfile[50];    /* Name for restart file */

    /* Names for the different uplift functions */
    String upliftName[] = {
        "plateau",
        "block",
        "tiltblock",
        "dynamic",
        "errfunc",
        "halfdome",
        "erfsub",
        "folds",
        "strike",
        "none",
    };

    /* Names for the different "east" boundary conditions */
    String eastBndCond[] = {
        "fixed zero",
        "fixed thickness",
        "no flux",
        "slope extrapolated",
        "none",
    };

    /* OPEN INPUT FILE */
    if( (ifp=fopen( infile, "r" ))==NULL )
        problem( "Unable to open input file." );

    /* READ PARAMETERS */
    end_time = ReadInt( ifp, "NT", 2 );
    Dt = ReadInt( ifp, "DT", 2 );
    dt_min = ReadDouble( ifp, "DTMIN", 5 );
    writeInterval = ReadInt( ifp, "WRTINT", 6 );
    dx = ReadDouble( ifp, "DX", 2 );

    optTransFn = ReadInt( ifp, "OPTTRANS", 5 );
    downstrmWt = ReadDouble( ifp, "DOWNWT", 6 );
    upstrmWt = 1.0 - downstrmWt;
    kf = ReadDouble( ifp, "KF", 2 );
    gkt = ReadDouble( ifp, "KT", 2 );
    mf = ReadDouble( ifp, "MF", 2 );
    nf = ReadDouble( ifp, "NF", 2 );
    gpf = ReadDouble( ifp, "PF", 2 );
    gtaucrit = ReadDouble( ifp, "TAUC", 4 );
    mb = ReadDouble( ifp, "MB", 2 );
    nb = ReadDouble( ifp, "NB", 2 );
    gkc = ReadDouble( ifp, "KC", 2 );   /* Channel width coef (10/97) */
    gkdp = ReadDouble( ifp, "KDP", 3 ); /* Channel depth coef (10/97) */
    gmdp = ReadDouble( ifp, "MDP", 3 ); /* Channel depth exponent (3/98) */

    optDiffuse = ReadInt( ifp, "OPTDIFF", 7 );
    kd = ReadDouble( ifp, "KD", 2 );

    mw = ReadDouble( ifp, "MW", 2 );   /* Weathering rate decay const */
    permChannelA = ReadInt( ifp, "CHANA", 5 );    /* Weathering threshold */
    precip_flow = ReadDouble( ifp, "PPT", 3 );    /* Mean precip intensity */
    gp_mean_ann = ReadDouble( ifp, "PMA", 3 );     /* Mean annual flow/area */
    optVarPrecip = ReadInt( ifp, "OPTP", 4 );     /* Precip variation option */
    if( optVarPrecip )
    {
        lambdap = ReadDouble( ifp, "LAMP", 4 );      /* Period of oscillation*/
        precip_var = ReadDouble( ifp, "PVAR", 4 );   /* Amount of variation */
    }

    optSat = ReadInt( ifp, "OPTSAT", 6 );
    if( optSat )
    {
        soilTrans = ReadDouble( ifp, "T", 1 );  /* Soil transmissivity */
        gwrchg = ReadDouble( ifp, "GWR", 3 );   /* Groundwater recharge rate */
    }

    solinit = ReadDouble( ifp, "C0", 2 );
    slopeCollapseActive = ReadInt( ifp, "OPTSLP", 6 );
    if( slopeCollapseActive ) 
    {
        soil_slope = ReadDouble( ifp, "SCR", 3 );
        optAreaDepLS = ReadInt( ifp, "OPTALS", 6 );     
        if( optAreaDepLS )
        {
            kls = ReadDouble( ifp, "KLS", 3 );
            angleOfRepose = ReadDouble( ifp, "ANGREP", 6 );
        }
        else angleOfRepose = soil_slope;
    }
    mode = ReadInt( ifp, "MODE", 4 );
    if( mode==RegionalScaleMode ) elevToTrack = ReadInt( ifp, "ZTR", 3 );
    else {
        kci = ReadDouble( ifp, "KCI", 3 );
        mci = ReadDouble( ifp, "MCI", 3 );
        nci = ReadDouble( ifp, "NCI", 3 );
        chanThreshold = ReadDouble( ifp, "TCI", 3 );
        chanThreshold = chanThreshold/kci;
    }
    chanarea = ReadDouble( ifp, "ACHAN", 4 );
    uptime = ReadInt( ifp, "UPDUR", 5 );
    if( uptime ) {
        upliftRate = ReadDouble( ifp, "UPRT", 4 );
        
        /* Parse the uplift type name and assign the proper code */
        fgets( tempstr, 80, ifp );
        fgets( tempstr, 80, ifp );
        uptype = 0;
        while( uptype < NumUpliftTypes && 
            !CompareStrings( tempstr, upliftName[uptype] ) ) uptype++; 

        /* Get parameters for specific uplift types */
        if( uptype==BLOCK || uptype==TILTBLOCK || uptype==FOLDS )
            fault = ReadInt( ifp, "FLOC", 4 );
        if( uptype==TILTBLOCK )
            fault_to_pivot = ReadDouble( ifp, "FTPIV", 5 );
        else if( uptype==BLOCK || uptype==FOLDS )
            vf = ReadDouble( ifp, "FPROP", 5 );
        if( uptype==ERRFUNC || uptype==ERFSUB || uptype==HALFDOME )
        {
            erfdist = ReadDouble( ifp, "EDIST", 5 );
            erfoffset = ReadInt( ifp, "EOFFS", 5 );
        }
        if( uptype==ERFSUB ) tsubstart = ReadInt( ifp, "SUBST", 5 );
        if( uptype==FOLDS ) {
            seedValue = ReadInt( ifp, "SEED", 4 );
            foldUplift( 0, INIT, seedValue, ifp );
        }
        if( uptype==STRIKESLIP )
        {
            fault = ReadInt( ifp, "FLOC", 4 );
            sliprate = ReadDouble( ifp, "SLIPRT", 6 );
        }
            
        /* Option for fault blocks along west boundary */
        ngrabens = ReadInt( ifp, "NGRAB", 5 );
        if( ngrabens )
        {
            grabenLength = ReadInt( ifp, "GRBLEN", 6 );
            if( grabenLength*ngrabens > XGRID )
                problem( "Too many grabens");
        }
    }

    /* Flexure option: if flexure is active, read flexural parameters */
    flexure = ReadInt( ifp, "OPTFLX", 6 );
    if( flexure )
    {
        flexr = ReadDouble( ifp, "FLEXR", 5 );
        rhoc = ReadDouble( ifp, "RHOC", 4 );
        rhom = ReadDouble( ifp, "RHOM", 4 );
        rhos = ReadDouble( ifp, "RHOS", 4 );
        basin_width = ReadDouble( ifp, "BASWID", 6 );
    }

    /* East boundary condition */
    fgets( tempstr, 80, ifp );
    fgets( tempstr, 80, ifp );
    east_boundary = 0;
    while( east_boundary < NumEastBdyOpts && 
        !CompareStrings( tempstr, eastBndCond[east_boundary] ) )
        east_boundary++;
    if( east_boundary==NumEastBdyOpts ) {
        printf( "East boundary type '%s' not recognized.\n", tempstr );
        printf( "Valid east boundary types are:\n" );
        for( i=0; i<NumEastBdyOpts; i++ )
            printf("\t%d\t%s\n", i, eastBndCond[i] );
        problem( "Invalid east boundary type." );
    } 

    /* Options for ponding, irregular boundary, discharge point sources,
     * and bedrock channel minimum gradient. */
    pond = ReadInt( ifp, "POND", 4 );
    optIrregularBoundary = ReadInt( ifp, "OPTBND", 6 );
    if( optIrregularBoundary ) 
    {
        fgets( tempstr, 80, ifp );
        fgets( boundaryfile, 50, ifp );
        boundaryfile[strlen(boundaryfile)-1] = '\0';
    }
    nPtSources = ReadInt( ifp, "NPT", 3 );
    if( nPtSources > 0 )
    {
        inputLoc = (int *)malloc( (unsigned) nPtSources*sizeof( int ) );
        Ainput = (int *)malloc( (unsigned) nPtSources*sizeof( float ) );
        for( i=0; i<nPtSources; i++ ) {
            inputLoc[i] = ReadInt( ifp, "ILOC", 4 );
            Ainput[i] = ReadInt( ifp, "AINP", 4 );
        }
    }
    optMinSlope = ReadInt( ifp, "OPTMS", 5 );
    if( optMinSlope ) minS = ReadDouble( ifp, "MINS", 4 );

    /* Names for auxiliary input files, output, & restart file (if any) */
    fgets( tempstr, 80, ifp );
    fgets( stratfile, 50, ifp );
    fgets( tempstr, 80, ifp );
    fgets( basefile, 50, ifp );
    fgets( tempstr, 80, ifp );
    fgets( outputfile, 50, ifp );
    restart = ReadInt( ifp, "OPTRST", 6 );
    if( restart )
    {
        fgets( tempstr, 80, ifp );
        fgets( restartfile, 50, ifp );
    }
    fclose( ifp );  
    /* Now we've finished reading the main input file! */

    /* APPEND NULL CHARACTER TO END OF FILENAMES */
    /* (NECESSARY WHEN STRINGS ARE READ FROM FILE) */
    stratfile[strlen(stratfile)-1] = '\0';
    outputfile[strlen(outputfile)-1] = '\0';
    basefile[strlen(basefile)-1] = '\0';
    restartfile[strlen(restartfile)-1] = '\0';

    /* Append code version and date stamp onto the input file. */
    ifp = fopen( infile, "a" );
    fprintf( ifp, "%s\n\t%s\n", VERSION, date() );
    fclose( ifp );

    /*  Set the transport function type. This is done by copying the 
     *  address of the appropriate function into the function pointer
     *  variable (*Transport)().
     */
    if( optTransFn==BAGNOLD )
        Transport = Bagnold;
    else if( optTransFn==THRESHOLD )
        Transport = ThresholdTransport;
    else if( optTransFn==MPM )
        Transport = MPMTransport;
    else
        Transport = GenericTransport;

    /*
     *  Initialize some variables.
     *      o cellarea: multiplying once now saves doing so repeatedly later.
     *      o chan: all cells initially are alluvial channels. 
     */
    cellarea = dx*dx;
    for( i=1; i<=XGRID; i++ )
        for( j=1; j<=YGRID; j++ )
            chan[i][j] = ALLUVIAL;

    /*
     *  Allocate space for arrays as needed.
     */
    if( slopeCollapseActive ) slp_fail=cmatrix(1,XGRID,1,YGRID);
    if( mode==RegionalScaleMode && elevToTrack==ChannelElev )
    {
        hillslopeReservoir=dmatrix(0,XGRIDP1,0,XGRIDP1); 
        meanElev=dmatrix(0,XGRIDP1,0,YGRIDP1);
    }
    if( optTransFn==BAGNOLD || optTransFn==THRESHOLD || optTransFn==MPM) 
        tau=dmatrix(0,XGRIDP1,0,YGRIDP1);

    /*
     *  Read in initial topography and stratigraphy
     *  If applicable, read boundary file.
     *      For now, the boundary is read as a binary file with the
     *  rows and cols reversed (the sort of file generated by PV-Wave).
     *  TODO: Make this cleaner, with a standard bin read method.
     */
    ReadLithology( stratfile, basefile, restartfile, restart );
    if( optIrregularBoundary )
    {
        if( (ifp=fopen(boundaryfile,"r"))==NULL ) {
            fprintf( stderr, "File: <%s>.\n", boundaryfile );
            problem( "Cannot open boundary data file.\n" );
        }
        fread( tmpCharArray, sizeof( tmpCharArray ), 1, ifp );
        fclose( ifp );
        for( i=0; i<=XGRIDP1; i++ ) 
            for( j=0; j<=YGRIDP1; j++ )
                boundary[i][j] = tmpCharArray[j][i];
    }
    
    /* If it's not a restart run, set the initial sediment thicknesses */
    if( !restart )
        for( i=1; i<=XGRID; i++ )
            for( j=1; j<=YGRID; j++ )
            {
                chansed[i][j] = solinit;
                sla[i][j] = solinit;
                slb[i][j] = solinit;
            }

    /* Here we establish all the cells along the west grid edge
     * as open boundaries. Cells along the east edge depend upon the
     * boundary condition: if it's a NO_FLUX boundary, they're closed
     * to sediment and water output; otherwise they're open.
     *  Also, if the grid is an irregular size, we set the elevation of
     * any open boundary cells now (as with any other cell, the elev is
     * the sum of the base elev plus sediment plus one or more rock layers).
     */
    for( j=0; j<=YGRIDP1; j++ )
        boundary[0][j] = OpenBoundary;
    if( east_boundary==NO_FLUX )
        for( j=0; j<=YGRIDP1; j++ ) boundary[XGRIDP1][j] = NoFluxBoundary;
    else
        for( j=0; j<=YGRIDP1; j++ ) boundary[XGRIDP1][j] = OpenBoundary;

    /* Now initialize the elevations of any open boundary cells within
     * the main part of the grid (this is only meaningful when an
     * irregular boundary is used). */
    for( i=1; i<=XGRID; i++ ) 
        for( j=1; j<=YGRID; j++ )
            if( boundary[i][j]==OpenBoundary )  
                elev[i][j] = 0;

    /* Adjustments for MPM transport function:
     * KF is assumed to be read in SI units, and is then modified to account
     * for:
     *  - rainfall intermittency factor
     *  - conversion from sed discharge to sed disch per unit cell area
     *  - inclusion of channel width constant terms: W = kc sqrt(Qm)
     *      (Qm = mean annual disch), = kc sqrt(Pm) sqrt(A)
     *      = kc sqrt(Pm) sqrt(cellarea) sqrt(N), N=area in cells (not m2)
     *  
*/
#define SECPERYEAR 31536000
    if( optTransFn==MPM )
    {
       /* Adjustments to KF: */
       kf *= gkc*sqrt( gp_mean_ann*cellarea ); /* channel width adjustment */
       kf *= gp_mean_ann / precip_flow; /* intermittency factor */
       kf = kf / cellarea;  /* converts to sed discharge per unit cell area */
       kf = SECPERYEAR * kf; /* converts from m/s to m/yr */
       /* Adjustments to KT: */
       gkt = (gkt*pow( gkc, -0.6 ) )  /* converts from shear per specific */
           / pow( gp_mean_ann, 0.3 ); /*  disch to shear per total Q */
       gkt = gkt * pow( precip_flow, 0.6 );  /* converts from shear per total*/
                                             /* disch to shear per area */
       precip_flow *= SECPERYEAR; /* convert precip to m/yr for bedrock erode*/
    }

    /* Precip flow rate in m3/yr: since precip_flow is read in as
     * a vertical accumulation rate (m/yr), we multiply by the cell
     * area to get the corresponding volumetric flow (which is the
     * discharge resulting from precipitation on one cell).
     * The variables precip_mean and precip_var are only used if
     * the option for time-varying precip (optVarPrecip) is active. 
     * Note: with MPM function, cell area is a separate term, not folded
     * into rainfall rate. Also, for MPM, units of precip are m/s NOT
     * m/yr! */
    precip_flow = precip_flow*cellarea;
    precip_mean = precip_flow;
    precip_var *= precip_mean;

    /*
     *  Scale the fluvial transport and erosion equations to convert
     *  sediment flux into elevation change. This is done by dividing 
     *  the erosion and transport coefficients by either the cell area
     *  or the channel area, depending on the simulation mode.
     *  (TODO: More careful parameterization of channel width needs
     *  to be worked out. For now, use channel area equal to cell
     *  to cell area, which is equivalent to assuming the whole cell
     *  lowers at the same rate as the channel.) 
     */
    if( mode==CatchmentMode || elevToTrack==AverageElev )
    {
        if( optTransFn!=MPM ) kf = kf/cellarea;
        for( k=1; k<=n_layers; k++ ) 
        {
            layer[k].Kmn *= (chanarea/cellarea);
            layer[k].K = layer[k].Kmn;
        }
    }
    else
        kf = kf/chanarea;
    kfmean = kf;          /* Mean transport coef, used only if precip varies */

    /* CALCULATE DIAGONAL DISTANCE BETWEEN GRID POINTS */
    dxx = sqrt( 2*dx*dx );

    /* CONVERT UPLIFT RATE FROM M/YR TO M/TIMESTEP */
    upliftRate *= Dt;

    /* 
     *  Set initial flexure values to zero. (This ensures that when 
     *  flexure is not active, elevations will not be affected by the 
     *  random contents of this array).
     */
    for( i=0; i<=XGRIDP1; i++ )
        w[i] = 0.0;

    /*
     *  Initialize some variables.
     */
    for( j=0; j<=YGRIDP1; j++ )
        elev[0][j] = 0.0;
    sed_yield = 0.0;

    /*
     *  Initialize tridiagonal matrices for ADI diffusion solver,
     *  if applicable.
     *  This is a time-for-space trade: we need to maintain NxMx6
     *  floating point arrays all the time, but don't need to reinitialize
     *  each time we call Hillslope.
     *  TODO: Fix ADI solver!
     *        Implement options for ADI vs explicit solvers
     */
/*  if( mode==CatchmentMode && kd>0.0 )
        ADIinit();
*/
    /* Set name for new restart file */
    strcpy( restart_name, outputfile );
    strcat( restart_name, ".restart" );

    /* Finally, call initialization functions to set up the list of cells,
     * open and initialize the output files, and if necessary, initialize
     * the flexure parameters.
     */
    MakeSortList();
    OpenOutputFiles( outputfile );
    if( flexure ) FlexInit( restart, solinit, basin_width );

} /**** end of initialize ****/



/*
 *  ADIinit:
 *  This function is out of commission until the ADI solver routine is
 *  repaired.
 */

void ADIinit()
{
    int k, r, c;

    /* Allocate memory space for the tridiagonal matrices:
     * Here uca, ucb and ucc are the upper, middle and lower
     * bands for one traverse, and v** the equivalent for the
     * other traverse.
     */
    uca = dvector(1,GRIDSIZE);
    ucb = dvector(1,GRIDSIZE);
    ucc = dvector(1,GRIDSIZE);
    vca = dvector(1,GRIDSIZE);
    vcb = dvector(1,GRIDSIZE);
    vcc = dvector(1,GRIDSIZE);

    /* Compute the R value, or more properly half the R value
     * since we're splitting the operation into two pieces
     * (one for each direction) every time step.
     */
    R = 0.5*kd*Dt/(dx*dx);

    /* Values for the coefficient matrices: upper and lower bands.
    for( k=1; k<=GRIDSIZE; k++ ) {
        uca[k] = -1.0;
        ucc[k] = -1.0;
        vca[k] = -1.0;
        vcc[k] = -1.0;
    }

    /* Values for the coefficient matrices: middle bands.
    for( k=1; k<=GRIDSIZE; k++ ) {
        ucb[k] = 1.0/R + 2.0;
        vcb[k] = 1.0/R + 2.0;
    }

    /* Replace these entries with 0 at the appropriate boundaries */
    for( r=1; r<=YGRID; r++ )   /* Rows for the odd matrix...*/
    {
        k = (r-1)*XGRID + 1;    /* Lefthand nodes */
        uca[k] = 0.0;
        k = r*XGRID;        /* Right hand nodes */
        ucc[k] = 0.0;
    }
    for( c=1; c<=XGRID; c++ )   /*...and columns for the even matrix */
    {
        k = (c-1)*YGRID + 1;    /* Top nodes */
        vca[k] = 0.0;
        k = c*YGRID;        /* Bottom nodes */
        vcc[k] = 0.0;
    }

}



/************************************************************************
 *                                                                      
 * LithologyVector:     Allocates a vector of Lithology structures      
 *                      with a range [nl..nh].                          
 *
 *
 *      Arguments: the range for the vector.
 *      Returns: a pointer to the allocated vector of 
 *               Lithology structures.
 *      Calls: (none)                                                   
 *      Called By: initialize                                           
 *                                                                      
/************************************************************************/

struct Lithology *LithVector(nl,nh)
int nl, nh;
{
    struct Lithology *v;

    v=(struct Lithology *)malloc((unsigned)(nh-nl+1)*sizeof(struct Lithology));
    if(!v) problem("allocation failure in LithVector()" );
    return v-nl;

}



/************************************************************************
 *
 *      ReadLithology:  This function reads the stratigraphy and base
 *                      files. If this is a restart from a previous 
 *                      run, the restart file is read instead of the
 *                      base file.
 *
 *      Arguments: the names of the stratigraphy, base topography, and
 *                 restart files, and a flag indicating whether or not
 *                 this is a restart of a previous run.
 *      Calls: problem, LithVector
 *      Called by: initialize
 *    
 ************************************************************************/

void ReadLithology( stratfile, basefile, restartfile, restart )
char *stratfile, *basefile, *restartfile;
int restart;
{
    int i, j, k,            /* counters for x & y coords */
            opt_vbl_thick;  /* Option for variable top layer thickness */
    FILE *tfp;              /* -> to initial topography file */
    float temp, temp2, temp3;

    /* Open the stratigraphy file */
    if( ( tfp=fopen( stratfile, "r" )) == NULL )
    {
        fprintf(stderr, "File <%s>\n", stratfile );
        problem( "Unable to open stratigraphy file." );
    }

    /* Read no. of layers, layered option, and allocate memory space
     *  for an array of Lithology records */
    fscanf( tfp, "%d", &n_layers );
    fscanf( tfp, "%d", &optLayeredLithology );
    layer = LithVector(1,n_layers);
    
    /* Here we read the variables for each lithology. 'K' is the bedrock
     * erodibility parameter. To avoid underflow, it is read in units
     * of itself x 10^-6. If a stream power erosion law is used, this
     * means it is read in units of mm-2 and converted to meters by
     * multiplying by 1,000,000. */
    for( k=1; k<=n_layers; k++ )
    {
        fscanf( tfp, "%f", &temp );
        layer[k].T = temp;
        fscanf( tfp, "%f", &temp );
        layer[k].Kmn = temp*0.0000010;
        layer[k].K = layer[k].Kmn; 
        fscanf( tfp, "%f", &temp );
        layer[k].W = temp;
        layer[k].W *= Dt;
        fscanf( tfp, "%f", &temp );
        layer[k].Scr = temp;
        fscanf( tfp, "%f", &temp );
        layer[k].tc = temp;
    }
    if( restart )
    {
        fclose( tfp );
        if( ( tfp=fopen( restartfile, "r" )) == NULL )
        {
            fprintf( stderr, "File <%s>\n", restartfile );
            problem( "Unable to open restart file." );
        }
        fscanf( tfp, "%d", &start_time );
        end_time += start_time;
        uptime += start_time;
        for( j=0; j<=YGRIDP1; j++ )
            for( i=0; i<=XGRIDP1; i++ )
            {
                fscanf( tfp, "%d%f%f%f", &rockType[i][j],&temp,&temp2,&temp3 );
                rock[i][j] = temp;
                base[i][j] = temp2;
                chansed[i][j] = temp3;
                sla[i][j] = temp3;
                slb[i][j] = temp3;
            }
        if( flexure )
        {
            fscanf( tfp, "%f", &temp );
            plat_thick = temp;
            fscanf( tfp, "%f", &temp );
            shelf_load = temp;
            for( i=0; i<=XGRIDP1; i++ )
            {
                fscanf( tfp, "%f", &temp );
                w_init[i] = temp;
                if( feof( tfp ) )
                    problem( "No flexure information in previous run.");
            }
        }

        /* If necessary, initialize the hillslopeReservoir array. 
           This variable records the volume of material above BASE
           contained in the hillslopes, which initially will be
           equal to the ROCK thickness plus any underlying layers,
           divided by the hillslope area (cell area - channel area). 
         */
        if( mode==RegionalScaleMode && elevToTrack==ChannelElev )
        {
            /*TODO: Open and read from hill reservoirs from file */
        }
    }
    else /* Not a restart run */
    {
        fscanf( tfp, "%d", &opt_vbl_thick );
        if( opt_vbl_thick )
            for( j=0; j<=YGRIDP1; j++ )
                for( i=0; i<=XGRIDP1; i++ )
                {
                    fscanf( tfp, "%d%f", &rockType[i][j], &temp );
                    rock[i][j] = temp;
                }
        else
            for( j=0; j<=YGRIDP1; j++ )
                for( i=0; i<=XGRIDP1; i++ )
                {
                    rockType[i][j]=1;
                    rock[i][j]=layer[1].T;
                }
        fclose( tfp );
        if( ( tfp=fopen( basefile, "r" )) == NULL )
        {
            printf( "File: <%s>\n",basefile );
            problem( "Unable to open base file." );
        }
        for( j=0; j<=YGRIDP1; j++ )
            for( i=0; i<=XGRIDP1; i++ )
            {
                fscanf( tfp, "%f", &temp );
                base[i][j]=temp;
            }

        /* If necessary, initialize the hillslopeReservoir array. 
         * This variable records the volume of material above BASE
         * contained in the hillslopes, which initially will be
         * equal to the ROCK thickness plus any underlying layers,
         * divided by the hillslope area (cell area - channel area). 
         */
        if( mode==RegionalScaleMode && elevToTrack==ChannelElev )
            for( j=0; j<=YGRIDP1; j++ )
                for( i=0; i<=XGRIDP1; i++ )
                {
                    hillslopeReservoir[i][j] = rock[i][j] + chansed[i][j];
                    if( optLayeredLithology )
                        for( k=rockType[i][j]; k<=n_layers; k++ )
                            hillslopeReservoir[i][j] += layer[k].T;
                    hillslopeReservoir[i][j] *= (cellarea-chanarea);
                } 
            
    }

    fclose( tfp );


} /**** end of ReadLithology ****/



/************************************************************************
 *                                                             
 *      MakeSortList:   This function creates a list of all non-boundary
 *                      grid cells. This list can then be sorted by
 *                      drainage area before computing fluvial erosion
 *                      and/or deposition in the Fluvial() function.
 *                      can be sorted according to drainage area.    
 *  
 ************************************************************************/

MakeSortList()
{
    int i, j, t=0;

    for( i=1; i<=XGRID; i++ )
        for( j=1; j<=YGRID; j++ ) 
            if( !boundary[i][j] )
            {
                cell[t].x = i;
                cell[t].y = j;
                t++;
            }
    numActiveCells = t;
    numChannelCells = numActiveCells;

} /* MakeSortList */



/************************************************************************
 *
 *      OpenOutputFiles:    This function opens the output files
 *                          and writes header information to them. 
 *
 *      Calls: (none)
 *      Called by:  initialize
 *
 ************************************************************************/

OpenOutputFiles( outfile )
char outfile[20];           /* Base name for output files, w/o ext */
{
    int ttime;              /* Total # of timesteps recorded */
    char fullname[20];      /* Full name of data files */

    ttime = (int) ( (float)(end_time-start_time) / (float)writeInterval );

    /* Current Time Step File */
    strcpy( timefile, outfile );
    strcat( timefile, ".time" );

    /* DRAINAGE NETWORK RECORDS */
    strcpy( fullname, outfile );
    strcat( fullname, ".net" );
    if( ( netfp=fopen(fullname,"wt") ) == NULL )
        problem("Cannot open output file");
    fprintf( netfp, "%d\n", ttime );

    /* ELEVATION RECORDS */
    strcpy( fullname, outfile );
    strcat( fullname, ".elev" );
    if( ( elevfp=fopen(fullname,"wt") ) == NULL )
        problem("Cannot open output file");
    fprintf( elevfp, "%d\n", XGRID+2 );
    fprintf( elevfp, "%d\n", YGRID+2 );
    fprintf( elevfp, "%f\n", dx );

    /* Mean ELEVATION RECORDS */
    if( mode==RegionalScaleMode && elevToTrack==ChannelElev ) {
        strcpy( fullname, outfile );
        strcat( fullname, ".mnelev" );
        if( ( mnelevfp=fopen(fullname,"wt") ) == NULL )
            problem("Cannot open output file");
        fprintf( mnelevfp, "%d\n", XGRID+2 );
        fprintf( mnelevfp, "%d\n", YGRID+2 );
        fprintf( mnelevfp, "%f\n", dx );
    }

    /* Sediment THICKNESS RECORDS */
    strcpy( fullname, outfile );
    strcat( fullname, ".sed" );
    if( ( solfp=fopen(fullname,"wt") ) == NULL )
        problem("Cannot open output file");

    /* DISCHARGE RECORDS */
    if( optSat )
      {
        strcpy( fullname, outfile );
        strcat( fullname, ".q" );
        if( ( qfp=fopen(fullname,"wt") ) == NULL )
            problem("Cannot open output file");
      }

    /* DRAINAGE AREA RECORDS */
    strcpy( fullname, outfile );
    strcat( fullname, ".area" );
    if( ( afp=fopen(fullname,"wt") ) == NULL )
        problem("Cannot open output file");

     /* SEDIMENT YIELD RECORDS */
    strcpy( fullname, outfile );
    strcat( fullname, ".sedyld" );
    if( ( sedfp=fopen(fullname,"wt") ) == NULL )
        problem("Cannot open output file");

    /* DENUDATION RATE RECORDS */
    strcpy( fullname, outfile );
    strcat( fullname, ".denud" );
    if( ( denudfp=fopen(fullname,"wt") ) == NULL )
        problem("Cannot open output file");

    /* SURFACE LITHOLOGY RECORDS */
    if( n_layers>1 ) {
        strcpy( fullname, outfile );
        strcat( fullname, ".sg" );
        if( ( lithfp=fopen(fullname,"wt") ) == NULL )
            problem("Cannot open output file");
        fprintf( lithfp, "%d\n", XGRID+2 );
        fprintf( lithfp, "%d\n", YGRID+2 );
        fprintf( lithfp, "%f\n", dx );
    }

    /* UPPERMOST LAYER ROCK THICKNESS RECORDS */
    if( n_layers>1 ) {
        strcpy( fullname, outfile );
        strcat( fullname, ".lyr1" );
        if( ( lyrfp=fopen(fullname,"wt") ) == NULL )
            problem("Cannot open output file");
        fprintf( lyrfp, "%d\n", XGRID+2 );
        fprintf( lyrfp, "%d\n", YGRID+2 );
        fprintf( lyrfp, "%f\n", dx );
    }

    /* CHANNEL TYPE RECORDS */
    strcpy( fullname, outfile );
    strcat( fullname, ".chan" );
    if( ( chanfp=fopen(fullname,"wt") ) == NULL )
        problem("Cannot open output file");

    /* SLOPE FAILURE RECORDS */
    if( slopeCollapseActive ) 
    {
        strcpy( fullname, outfile );
        strcat( fullname, ".slp" );
        if( ( slpfp=fopen(fullname,"wt") ) == NULL )
            problem("Cannot open output file");
    }

    /* FLEXURE RECORDS */
    if( flexure )
    {
        strcpy( fullname, outfile );
        strcat( fullname, ".flex" );
        if( ( flexfp=fopen(fullname,"wt") ) == NULL )
            problem("Cannot open output file");
    }

    /* SHEAR STRESS RECORDS */
    if( optTransFn==BAGNOLD || optTransFn==THRESHOLD || optTransFn==MPM )
    {
        strcpy( fullname, outfile );
        strcat( fullname, ".tau" );
        if( ( taufp=fopen(fullname,"wt") ) == NULL )
            problem("Cannot open output file");
    }

    /* Write header information to the relevant output files */
    fprintf( solfp, "%d\n", XGRID+2 );
    fprintf( solfp, "%d\n", YGRID+2 );
    fprintf( solfp, "%f\n", dx );
    fprintf( afp, "%d\n", XGRID+2 );
    fprintf( afp, "%d\n", YGRID+2 );
    fprintf( afp, "%f\n", dx );
    if( optSat )
      {
        fprintf( qfp, "%d\n", XGRID+2 );
        fprintf( qfp, "%d\n", YGRID+2 );
        fprintf( qfp, "%f\n", dx );
      }
    fprintf( denudfp, "%d\n", XGRID+2 );
    fprintf( denudfp, "%d\n", YGRID+2 );
    fprintf( denudfp, "%f\n", dx );
    fprintf( chanfp, "%d\n", XGRID+2 );
    fprintf( chanfp, "%d\n", YGRID+2 );
    fprintf( chanfp, "%f\n", dx );
    if( optTransFn==BAGNOLD || optTransFn==THRESHOLD || optTransFn==MPM )
    {
        fprintf( taufp, "%d\n", XGRID+2 );
        fprintf( taufp, "%d\n", YGRID+2 );
        fprintf( taufp, "%f\n", dx );
    }

} /**** end of OpenOutputFiles ****/



/************************************************************************
 *
 *  FlexInit:   This function initializes the flexural isostasy
 *              equations. These equations assume an elastic
 *              beam of uniform rigidity, oriented along the model
 *              X axis (that is, perpendicular to the lefthand
 *              boundary). For a discussion of the equations see
 *              Turcotte and Schubert (1983). This function assumes
 *              that the initial crustal load is in isostatic 
 *              equilibrium; this is done by calculating the initial
 *              flexural response and then subtracting this from
 *              the initial and all subsequent solutions.
 *
 *      Calls:  (none)
 *      Called by:  initialize
 *
 ************************************************************************/

void FlexInit( restart, solinit, l )
int restart;    /* flag indicating restart of an existing run */
float solinit;  /* initial soil thickness */
double l;       /* initial basin width (region over which seds spread) */
{
    int i, j;       /* Counters */
    double x,       /* Distance from a loading point (m) */
        w0,     /* Deflection (m) at a loading point */
        wn;     /* Deflection at a distance x from a loading pt */

    /* Initialize variables */
    for( i=0; i<=XGRIDP1; i++ ) {
        w[i] =0.0;
        V[i] = 0.0;
        plat_w[i] = 0.0;
        shelf_w[i] = 0.0;
    }

    /*  Calculate flexural parameter ALPHA. Also calculate UNITLOAD
     *  as the flexural response to unit load (not including gravity
     *  or density or width-of-load terms).
     */
    alpha = (4*flexr)/(rhom*GRAV);
    alpha = pow(alpha,0.25) * 10000.0;
    unitload = (pow(alpha,3.0)*0.00000001)/(8.0*flexr);

    /*
     *  Add a "fictitious" extra 500km for flexure of continuous
     *  plateau to the west....
     */
    if( !restart )
    {
        shelf_load = 0.0;
        plat_thick = 0.0;       
        for( i=1; i<=n_layers; i++ )
            plat_thick += layer[i].T;
        plat_thick += solinit*(rhos/rhoc);
    }
    w0 = unitload*rhoc*GRAV*5000.0*0.00000001;

        for( i=1; i<=100; i++ )
        {
                for( j=1; j<=XGRIDP1; j++ ) {
                        x = 2500.0+(i-1)*5000.0 + dx*(0.5+(XGRID-j));
                        plat_w[j] += w0*exp(-x/alpha)*(cos(x/alpha)+sin(x/alpha)
);
                }
        }

        /*
         *  Parameter for converting sediment flux to offshore sediment
         *  pile. Pile is L meters long.
         */
        shelf_param = ((dx*rhos*GRAV*5000.0)/(l*(double)YGRID));

        /*
         *  Compute flexural response for unit load on shelf.
         */
        for( i=1; i<=200; i++ )
        {
                for( j=0; j<=XGRIDP1; j++ ) {
                        x = i*5000.0 + dx*j;
                        shelf_w[j] += unitload*exp(-x/alpha)*(cos(x/alpha)+sin(x
/alpha));
                }
        }

        /*
         *  Calculate an initial deflection. This profile will be subtracted
         *  from the new profile for each time step, reflecting the assumption
         *  that the initial mass distribution is in equilibrium.
         *      First, the deflection due to the load in the area modeled is
         *  calculated. Then, the deflection due to the load of the interior
         *  plateau (=response to unit plateau thickness x plateau thickness)
         *  is added to this profile.
         */
        if( !restart )
        {
                for( i=0; i<=XGRIDP1; i++ ) w_init[i] = 0.0;
                w0 = unitload*plat_thick*rhoc*GRAV*dx*0.00000001;
                for( i=1; i<=XGRIDP1; i++ ) {
                        for( j=0; j<=XGRIDP1; j++ ) {
                                x = abs(j-i) * dx;
                                wn = w0*exp(-x/alpha)*(cos(x/alpha)+sin(x/alpha)
);
                                w_init[j] += wn;
                        }
                }
        }
        for( i=0; i<=XGRIDP1; i++ )
                plat_init[i] = plat_w[i]*plat_thick;

} /**** end of FlexInit ****/



/************************************************************************
 *
 *  ResetVariables: This function is called at the beginning of each
 *                  time step to recompute cell elevations as the sum
 *                  of the base elevation plus the overlying strata
 *                  and sediment, plus any flexural adjustments.
 *                  This is where uplift, which modifies the basal
 *                  surface, gets translated into elevation. This
 *                  function also adjusts the model boundaries for
 *                  the new time step (either directly or by calling
 *                  one or more helper functions).
 *
 *      Calls:  adjustEastBoundary, wrapNSBoundaries
 *      Called by:  main
 *
 ************************************************************************/

void ResetVariables()
{
    int i, j, k;                    /* counters */
    double sum,                     /* accum'd sum of bdy elev'n */
        cellareaInverse = 1.0/cellarea; /* 1/(area of a cell in sq m) */

    /*
     *  Compute elevations for the new time step as a function of
     *  total rock layer thickness and soil thickness, plus modification
     *  due to lithosphere flexure if appropriate. Also reset flood
     *  status flags and discharges.
     *      And while we're at it, compute mean denudation.
     */
    denud = 0.0;
    for( i=1; i<=XGRID; i++ )
        for( j=1; j<=YGRID; j++ )
        {
           denud += delz[i][j][0];
            if( !boundary[i][j] )
            {
                elev[i][j] = base[i][j] + rock[i][j] + chansed[i][j] - w[i];
                if( optLayeredLithology )
                    for( k=rockType[i][j]+1; k<=n_layers; k++ )
                        elev[i][j] += layer[k].T;
                flood[i][j] = FALSE;
            } /* for j */
        }
    denud = denud/(float)GRIDSIZE;

    /*
     *  If we're tracking the channel bed elevations as ELEV in
     *  RegionalScaleMode, then we'll also want to know the mean cell
     *  elevations. We calculate that here as the total vol of material
     *  in the cell divided by its area. Volume is hillslopeReservoir
     *  plus channel material (rock + chansed) times channel area.
     */ 
    if( mode==RegionalScaleMode && elevToTrack==ChannelElev )
        for( i=1; i<=XGRID; i++ )
            for( j=1; j<=YGRID; j++ )
                meanElev[i][j] = cellareaInverse*(hillslopeReservoir[i][j]+
                                 chanarea*(rock[i][j]+chansed[i][j])) - w[i];

    /*
     *  Set east (righthand) boundary elevation
     */
    if( east_boundary==FIXED_ZERO )
        for( j=0; j<=YGRIDP1; j++ ) elev[XGRIDP1][j] = 0;
    else if( east_boundary==FIXED_THICKNESS )
        for( j=0; j<=YGRIDP1; j++ )
        {
            elev[XGRIDP1][j] = base[XGRIDP1][j] + rock[XGRIDP1][j] + chansed[XGRIDP1][j] - w[i];
            if( optLayeredLithology )
                for( k=rockType[i][j]+1; k<=n_layers; k++ )
                    elev[XGRIDP1][j] += layer[k].T;
        }
    else if( east_boundary==NO_FLUX )
        for( j=0; j<=YGRIDP1; j++ ) elev[XGRIDP1][j] = elev[XGRID][j];
    else adjustEastBoundary();

    /*
     *      With half dome uplift, set floating west boundary.
     */
    if( uptype==HALFDOME )
    {
        for( j=1; j<=YGRID; j++ )
            sum += elev[1][j];
        sum = sum/(float)YGRID;
        for( j=0; j<=YGRIDP1; j++ )
            elev[0][j] = sum;
    }

    for( i=0; i<=XGRIDP1; i++ )
        for( j=0; j<=YGRIDP1; j++ )
            load[i][j] = 0.0;

    /* WRAP NORTH & SOUTH BOUNDARIES */
    wrapNSBoundaries();

} /**** end of ResetVariables ****/



/************************************************************************
 *
 *  wrapNSBoundaries:   This function sets the elevations and other
 *                      properties along the north and south (top 
 *                      and bottom) boundaries. To make these boundaries
 *                      reflective, the values along the bottom 
 *                      boundary are set equal to the values along the
 *                      uppermost (non-boundary) grid row, and vice-
 *                      versa.
 *
 *      Calls:  (none)
 *      Called by:  ResetVariables
 *
 ************************************************************************/

void wrapNSBoundaries()
{
    int i;

    for( i=1; i<=XGRID; i++ )
    {
        elev[i][0] = elev[i][YGRID];
        chansed[i][0] = chansed[i][YGRID];
        base[i][0] = base[i][YGRID];
        rock[i][0] = rock[i][YGRID];
        elev[i][YGRIDP1] = elev[i][1];
        chansed[i][YGRIDP1] = chansed[i][1];
        base[i][YGRIDP1] = base[i][1];
        rock[i][YGRIDP1] = rock[i][1];
    }

}



/************************************************************************
 *
 *  ShiftAverages:  This function maintains a three-point moving
 *                  average for sediment thickness and denudation
 *                  rate. It simply shifts values for each of the
 *                  past 3 time steps downward by one, making room
 *                  for the new time step. This averaging is only
 *                  for reporting purposes and does not affect the
 *                  model's computations.
 *
 *      Calls:  (none)
 *      Called by:  main
 *
 ************************************************************************/

void ShiftAverages()
{
    int i, j;
        
    sed_yield = 0.0;
    for( i=1; i<=XGRID; i++ )
        for( j=1; j<=YGRID; j++ )
        {
            delz[i][j][2] = delz[i][j][1];
            delz[i][j][1] = delz[i][j][0];
            delz[i][j][0] = 0.0;
            slb[i][j] = sla[i][j];
            sla[i][j] = (float)chansed[i][j];
        }

}



/************************************************************************
 *
 *  adjustEastBoundary: This function is only called when the user
 *                      wants to implement a special kind of righthand
 *                      (east) boundary that allows for a finite surface
 *                      slope AWAY from the grid, but does not fix
 *                      the boundary elevation. The function finds
 *                      the average gradient along the 2 grid columns
 *                      adjacent to the righthand boundary, and then
 *                      sets the elevation along this boundary such
 *                      that the AVERAGE from edge to boundary is 
 *                      equal to the average immediately inboard.
 *                      This technique is useful in cases where the
 *                      landscape tilts to the right at a wavelength
 *                      larger than the width of the model grid.
 *
 *      Calls:  (none)
 *      Called by:  ResetVariables
 *
 ************************************************************************/

void adjustEastBoundary()
{
    int j;
    double z_avg=0.0, dz_avg=0.0, z0;

    for( j=1; j<=YGRID; j++ )
    {
        z_avg += elev[XGRID][j];
        dz_avg += (elev[XGRID][j]-elev[XGRID-2][j]);
    }
    z0 = (z_avg+0.5*dz_avg)/(double)YGRID;
    for( j=0; j<=YGRIDP1; j++ )
        elev[XGRIDP1][j] = z0;

}



/************************************************************************
 *  
 *      FindNextCell:   Calculates the drainage direction for cell (x,y).
 *                      The function does this by computing the slope
 *                      in each direction, and identifying the steepest
 *                      of these. The adjacent cell lying in that
 *                      direction is recorded as the neighbor cell of
 *                      (x,y). The algorithm checks to see whether a
 *                      change of flow direction has occurred; if it
 *                      has, the flag streamCaptureOccurred is set to
 *                      TRUE, indicating that flow paths will need to
 *                      be recomputed. The algorithm also checks 
 *                      whether the cell forms a depression. If so,
 *                      the cell's flood status is set to BASIN.
 *
 *      Calls:  (none)
 *      Called by:  main
 *
 ************************************************************************/

void FindNextCell( x, y )
int x;
int y;
{
    int i,j, yy,
         ox, oy;    /* Previous drainage dir */
    double odx, os,
        slp[3][3],
        S;              /* Steepest slope at current cell       */

    /* Store the previous drainage direction and slope distance */
    ox = nbrx[x][y];
    oy = nbry[x][y];
    odx = delx[x][y];

    /* To start, assume that the cell is not a closed basin */
    flood[x][y] = FALSE;

    /* Compute slope to each surrounding cell. For diagonals,
     * multiply by SQRT(2).
     * TODO: Add alternative fns for regional mode */
    for( j= -1; j<=1; j++ )
    {
        yy=y+j;
        if( yy==YGRIDP1 ) yy=1;
        else if( yy==0 ) yy=YGRID;
        for( i= -1; i<=1; i++ )
            slp[i+1][j+1] = elev[x][y]-elev[x+i][yy];
    }
    slp[2][2] *=0.7071;
    slp[2][0] *=0.7071;
    slp[0][2] *=0.7071;
    slp[0][0] *=0.7071;

    /* Begin with null values */ 
    S = -100.0;
    nbrx[x][y] = 0;
    nbry[x][y] = 0;

    /* Find the greatest slope and assign the neighbor cell.
     * We only consider cells that aren't no-flux boundary cells.
     */
    for( i= -1; i<=1; i++ )
        for( j= -1; j<=1; j++ )
        {
            if( slp[i+1][j+1] > S && boundary[x+i][y+j]!=NoFluxBoundary )
            {
                S=slp[i+1][j+1];
                nbrx[x][y] = x+i;
                nbry[x][y] = y+j;
                delx[x][y]=(i==0 || j==0) ? dx : dxx;
            }
        }

    if( nbry[x][y]==YGRIDP1 ) nbry[x][y]=1;
    else if( nbry[x][y]==0 ) nbry[x][y]=YGRID;

    if( ( ox != nbrx[x][y] ) || ( oy != nbry[x][y] ) )
        streamCaptureOccurred = TRUE;

    /* IF STEEPEST SLOPE IS <= 0 THEN CELL IS A CLOSED BASIN */
    if( S <= 0 )
    {
        delx[x][y] = dx;
        slope[x][y] = 0;
        nbrx[x][y] = -1;
        nbry[x][y] = -1;
        flood[x][y] = BASIN;
    }

    /* Otherwise, assign greatest slope to current cell's 'slope' field */
    slope[x][y] = S/dx;

    /*Xif( delx[x][y]==0 && S > 0 ) 
        printf( "DELX is ZERO!!\n");*/

} /* FindNextCell */



/************************************************************************
 *
 *  DoPonds:    This function implements the algorithm that handles
 *              the drainage of sinks, that is, cells that have no
 *              downslope drainage. The algorithm works by finding the
 *              lowest cell along the perimeter of a sink cell or a
 *              group of sink cells. This cell is tested for downslope
 *              drainage in any direction other than back toward the
 *              flooded cells. If this cell has drainage, it is marked
 *              as the exit point for the sink cell (or group of sink
 *              cells). Otherwise, the cell is added to the flooded
 *              area and the process is repeated.
 *                  Because regions of closed drainage can consist of
 *              many cells, a linked list is used to store the grid
 *              coordinates of each cell in the flooded region. The
 *              "flood status" is encoded in the global grid variable
 *              FLOOD. A cell that is not UNFLOODED may be any of the
 *              following:
 *
 *                  BASIN---denotes a sink cell whose drainage has not
 *                          been determined. If the ponding algorithm
 *                          is not invoked, all BASIN cells are treated
 *                          as "black holes" in which all water evaporates
 *                          and all sediment is trapped. The ponding
 *                          algorithm converts all BASINs to OLDPONDs.
 *                  NEWPOND-denotes a cell that is part of a flooded
 *                          region, and whose ultimate drainage has
 *                          not yet been determined (ie, it is still
 *                          being processed). This distinguishes
 *                          cells that are part of the current flooded
 *                          region from areas whose drainage has
 *                          already been found. 
 *                  OLDPOND-In flat, irregular landscapes, different
 *                          flooded areas will often coalesce. OLDPOND
 *                          denotes a cell whose drainage has already
 *                          been computed during a previous pass through
 *                          the algorithm. When the algorithm encounters
 *                          an OLDPOND cell, it will add it to the group
 *                          of currently flooded cells and re-code it
 *                          as a NEWPOND. When drainage has been found
 *                          for every sink, all flooded cells will be
 *                          coded as OLDPONDs.
 *
 *      Calls:  Addlist, FindDryNeighbor
 *      Called by:  main
 *
 ************************************************************************/

DoPonds( cnx, cny )
int cnx, cny;
{
    struct PondCell *top, *bottom, *curr;
    int n, a, b, bb, ncells;
    struct CellPoint lowcell;
    char done = FALSE;

    /*printf("DoPonds()...");*/

    /* Hack: don't bother with boundary nodes. TODO: fix this properly */
    if( optIrregularBoundary )
        if( boundary[cnx][cny] ) return;

    /* CREATE FIRST ELEMENT IN LIST OF FLOODED CELLS */
    top = AddList();

    /* current cell is flooded */
    flood[cnx][cny] = NEWPOND;
    top->x = cnx;
    top->y = cny;
    top->next = NULL;
    bottom = top;
    ncells = 1;

    while( !done )
    {
        curr = top;             /* start @ top of list */
        lowcell.z = 1e9;      /* everything's lower than this */

        /* find the lowest cell on the perimeter of the pond */
        for( n=1; n<=ncells; n++ )
        {
            for( a=curr->x - 1; a<=curr->x+1; a++ )
                for( bb=curr->y - 1; bb<=curr->y+1; bb++ )
                {
                    /* IF N/S BOUNDARY, WRAP AROUND */
                    if( bb==0 ) b = YGRID;
                    else if( bb==YGRIDP1 ) b = 1;
                    else b = bb;
                                        
                    /* if cell isn't already flooded, test its elevation */
                    if( (flood[a][b]==FALSE || flood[a][b]==BASIN) 
                            && boundary[a][b]!=NoFluxBoundary )
                    {
                        if( elev[a][b] < lowcell.z )
                        {
                            lowcell.x = a;
                            lowcell.y = b;
                            lowcell.z = elev[a][b];
                        }
                    }

                    /* if cell is flooded but not on list, add it to list */
                    else if( flood[a][b] == OLDPOND )
                    {
                        bottom->next = AddList();
                        bottom = bottom->next;
                        bottom->x = a;
                        bottom->y = b;
                        bottom->next = NULL;
                        flood[a][b] = NEWPOND;
                        ncells++;
                    }
                } /* for b */
            curr = curr->next;      /* next cell on list */
        }  /* for ncells */

        /* if lowcell is an outflux boundary cell, call it an outlet */
        /* (Note: if through input error the elevations are so high */
        /* that no outlet is found, model may crash on this line    */
        /* bcs lowcell.x, .y is uninitialized)                      */
        if( boundary[lowcell.x][lowcell.y]==OpenBoundary ) done = TRUE;

        /* if not, find lowest unflooded cell surrounding LOWCELL */
        else
        {
            /* If LOWCELL has drainage, LOWCELL is pond's outlet */
            if( FindDryNeighbor( lowcell.x, lowcell.y ) ) 
                done = TRUE;
                        
            /* ...OTHERWISE, ADD LOWNODE TO THE LIST AND CONTINUE */
            else
            {
                bottom->next = AddList();     /* add new cell */
                bottom = bottom->next;
                bottom->x = lowcell.x;        /* give it lowcell coords. */
                bottom->y = lowcell.y;
                bottom->next = NULL;
                flood[lowcell.x][lowcell.y] = NEWPOND;   /* lowcell flooded */
                ncells++;
            }

        } /* else */

    } /* while */
        
    /* Now we've found an outlet: flood all cells and erase list */
    curr = top;
    for( n=1; n<=ncells; n++ )
    {
        nbrx[curr->x][curr->y] = lowcell.x;
        nbry[curr->x][curr->y] = lowcell.y;
        flood[curr->x][curr->y] = OLDPOND;
        top = curr;
        curr = curr->next;
        free( top );
    }

    /*printf("done\n");*/
    

} /**** end of DoPonds ****/



/************************************************************************
 *  
 *  AddList:    This function allocates memory for a new PondCell 
 *              structure and returns a pointer to the new structure.
 *              It is used to create new list elements in the linked
 *              list of flooded cells used by DoPonds.
 *
 *      Calls:  problem
 *      Called by:  DoPonds
 *
 ************************************************************************/

struct PondCell *AddList()
{
    struct PondCell *newCell;

    newCell = (struct PondCell *)malloc( (unsigned)sizeof( struct PondCell ) );
    if( newCell==NULL )
        problem( "Can't allocate memory for cell in ADDLIST\n" );

    return( newCell );
}


/************************************************************************
 *  
 *  FindDryNeighbor:    This function is used to find an outlet for
 *                      a cell on the periphery of a flooded region.
 *                      It's identical to FindNextCell, except that
 *                      in order to be considered as a valid outlet,
 *                      a potential flow direction must lead to a
 *                      cell that isn't itself already part of a
 *                      flooded region. The return value of the
 *                      function indicates whether or not it has
 *                      succeeded in finding a flow direction that
 *                      leads downhill to an unflooded ("dry")
 *                      neighbor.
 *
 *      Calls:  isdry
 *      Called by:  DoPonds
 *
 ************************************************************************/

int FindDryNeighbor( x, y )
int x;
int y;
{
    int i,j, yy,
         ox, oy;    /* Previous drainage dir */
    double odx, os,
        slp[3][3],
        S;              /* Steepest slope at current cell       */

    ox = nbrx[x][y];
    oy = nbry[x][y];
    odx = delx[x][y];

    /* Compute slope to each surrounding cell */
    /* TODO: Add alternative fns for regional mode */
    for( j= -1; j<=1; j++ )
    {
        yy=y+j;
        if( yy==YGRIDP1 ) yy=1;
        else if( yy==0 ) yy=YGRID;
        for( i= -1; i<=1; i++ )
            slp[i+1][j+1] = elev[x][y]-elev[x+i][yy];
    }
    slp[2][2] *=0.7071;
    slp[2][0] *=0.7071;
    slp[0][2] *=0.7071;
    slp[0][0] *=0.7071;

    /* Begin with null values */ 
    S = 0.0;
    nbrx[x][y] = 0;
    nbry[x][y] = 0;

    /* Find the greatest slope and assign the neighbor cell.
     * We only consider cells that aren't no-flux boundary cells.
     * Only care about dry cells!
     */
    for( i= -1; i<=1; i++ )
        for( j= -1; j<=1; j++ )
        {
            yy=y+j;
            if( yy==YGRIDP1 ) yy=1;
            else if( yy==0 ) yy=YGRID;
            if( slp[i+1][j+1] > S && boundary[x+i][yy]!=NoFluxBoundary && isdry( x+i, yy, x, y ) )
            {
                S=slp[i+1][j+1];
                nbrx[x][y] = x+i;
                nbry[x][y] = yy;
                delx[x][y]=(i==0 || j==0) ? dx : dxx;
            }
        }

    /* Otherwise, assign greatest slope to current cell's 'slope' field */
    if( S <= 0.0 ) return( 0 );
    slope[x][y] = S/delx[x][y];
    return( 1 );

} /* FindNextCell */



/************************************************************************
 *  
 *  isdry: This function reports whether a given cell (x,y) is flooded.
 *         If the flood status code is zero, or BASIN, the cell either
 *         has drainage or has not yet been processed by the ponding
 *         algorithm, in which case it's "dry."  
 *         
 *      Calls:  (none) 
 *      Called by:  FindDryNeighbor 
 *
 ************************************************************************/

int isdry( x, y, lowx, lowy )
int x, y, lowx, lowy;
{
    int outx, outy;
        
    if( !flood[x][y] || flood[x][y]==BASIN )
        return( TRUE );

    if( flood[x][y] == NEWPOND )
        return( FALSE );

    /* if it passes these tests, it's an old pond                     */
    /* occasionally, an "old pond" cell can have an outlet that whose */
    /* elevation is exactly equal to the current "lowcell's" elev     */
    /* (or it may even be the "lowcell")                              */
    /* In this situation, the test fails, since it would allow for    */
    /* the possibility of endless loops.                              */
    outx = nbrx[x][y];
    outy = nbry[x][y];
    if( elev[outx][outy]==elev[lowx][lowy] )
        return( FALSE );
    else
        return( TRUE );

}



/************************************************************************
 *    
 *      StreamTrace:    Computes area drained/total discharge by
 *                      tracing the drainage route of precipitation 
 *                      entering each cell.
 *
 *        Calls: (none)
 *        Called by: main
 *        Modifications: 3/98 added loop to compute channel depth
 *
 ************************************************************************/

StreamTrace()
{
    int i, j,   /* counters for x & y coords    */
        p, q,   /* x & y coords of cells traced along stream line */ 
        newp,   /* used so that old value of p won't get lost   */
        test;   /* counter tests for endless loops when cells   */
                /*      point to one another.                   */

    for( i=0; i<=XGRIDP1; i++ )
        for( j=0; j<=YGRIDP1; j++ ) 
            area[i][j] = 0;

    for( i=1; i<=XGRID; i++ )
        for( j=1; j<=YGRID; j++ ) 
        if( !boundary[i][j] )
        {
            p = i;
            q = j;
            area[i][j]++;
            test = 0;
            while( boundary[p][q]!=OpenBoundary && flood[p][q]!=BASIN )
            {
                test++;
                if( test > 1000 ) bugdump( p, q );
                newp = nbrx[p][q];
                q = nbry[p][q];
                p = newp;
                area[p][q]++;
            } /* while */
        } /* for j */

    /* Add the effect of any discharge point sources along the east bdy */
    for( i=0; i<nPtSources; i++ )
    {
        p = XGRID;
        q = inputLoc[i];
        area[p][q] += Ainput[i];
        while( boundary[p][q]!=OpenBoundary && flood[p][q]!=BASIN )
        {
            newp = nbrx[p][q];
            q = nbry[p][q];
            p = newp;
            area[p][q] += Ainput[i];
        }
    }

    if( gkdp>0 )
        for( i=1; i<=XGRID; i++ )
            for( j=1; j<=YGRID; j++ )
                chanDepth[i][j] = gkdp*pow( (double)area[i][j], gmdp );

} /**** end of StreamTrace ****/



/************************************************************************
 *    
 *  CalcDischarge
 *
 *  Computes discharge as a function of drainage area and, if the
 *  saturation option is "on", the groundwater conveyance capacity.
 *
 *        Calls: (none)
 *        Called by: main
 *        Created: 6/97
 *
 ************************************************************************/

CalcDischarge()
{
    int i, j,   /* counters for x & y coords    */
        p, q,   /* x & y coords of cells traced along stream line */ 
        newp,   /* used so that old value of p won't get lost   */
        test;   /* counter tests for endless loops when cells   */
                /*      point to one another.                   */
    double excessCap,
           locRunoff;

    if( !optSat )
      {
        for( i=0; i<=XGRIDP1; i++ )
            for( j=0; j<=YGRIDP1; j++ )
                disch[i][j] = area[i][j]*precip_flow;
      }
    else
      {
        for( i=1; i<=XGRID; i++ )
            for( j=1; j<=YGRID; j++ ) 
                disch[i][j] = 0;
         for( i=1; i<=XGRID; i++ )
            for( j=1; j<=YGRID; j++ ) 
                if( !boundary[i][j] )
                {
                    excessCap = soilTrans*slope[i][j] - gwrchg*area[i][j];
                    if( excessCap<0 ) excessCap = 0;
                    locRunoff = precip_flow - excessCap;
		    /*printf("(%d %d) ec: %f  lr: %f\n",i,j,excessCap,locRunoff);*/
                    if( locRunoff > 0 ) /* There is runoff; route downstrm*/
		      {
                        disch[i][j] += locRunoff;
                        p = i;
                        q = j;
                        while( boundary[p][q]!=OpenBoundary && flood[p][q]!=BASIN )
                        {
                            newp = nbrx[p][q];
                            q = nbry[p][q];
                            p = newp;
                            disch[p][q]+=locRunoff;
		        }
		      }
		} /* if */
        } /* else */
    /*         for( i=1; i<=XGRID; i++ )
            for( j=1; j<=YGRID; j++ ) 
                if( !boundary[i][j] && area[i][j]>1 )
printf("(%d %d) a: %d  q: %f\n",i,j,area[i][j],disch[i][j]/10000);*/
}



/************************************************************************
 *
 *  SedFluxToChannels:  This function is meant to be developed as part
 *                      of a parameterization scheme for sub-grid-scale
 *                      sediment delivery in RegionalScaleMode.
 *
 ************************************************************************/

SedFluxToChannels()
{
/* TODO: Add some function here */
}



/************************************************************************
 *      Function WEATHER
 ************************************************************************
 *    
 *  Weather:    This function calculates the addition to soil thickness
 *              during a time step, using an analytical solution to the
 *              equation  dC/dt = W * exp( -C ). With large soil thick-
 *              nesses, numerical errors in computing e raised to a
 *              large number can actually result in decreases in the
 *              computed soil thickness (and increases in the rock
 *              layer thickness). To avoid this error, the
 *              computation is only performed for soils less than 100m
 *              thick.
 *                  Modified for version 5.10 to allow weathering only
 *              where drainage area, in cells, is less than a specified
 *              minimum threshold.
 *
 *      Calls: (none)
 *      Called by: main (in catchment-scale mode only)
 *
 ************************************************************************/

Weather()
{
    int i, j; 
    double oldsed,
        kwdivmw;

    for( i=1; i<=XGRID; i++ )
        for( j=1; j<=YGRID; j++ )
        {
            if( area[i][j]<permChannelA && chansed[i][j] <= 100.0 ) {
                kwdivmw = layer[rockType[i][j]].W/mw;
                oldsed = chansed[i][j];
                chansed[i][j] = mw*log( exp(chansed[i][j]/mw) + kwdivmw );
                rock[i][j] -= (chansed[i][j]-oldsed);
                if( rock[i][j] < 0.0 )
                {
                    rockType[i][j]++;
                    rock[i][j] += layer[rockType[i][j]].T;
                }
            }
        }
}



/************************************************************************
 *    
 *  FindChannels:   Called only in catchment-scale mode, and then only
 *                  when a channelization threshold is used, this 
 *                  function tests each cell for the presence of a
 *                  stream channel by comparing the value of a function
 *                  kc Q^mci S^nci with a specified threshold value.
 *                  To speed up computation, the algorithm first
 *                  tests for the special case in which nci=0 and 
 *                  mci=1.
 *
 *      Calls: (none)
 *      Called by: main (in catchment-scale mode only)
 *
 ************************************************************************/

FindChannels()
{
    int i, j;
    double fn, S;

    numChannelCells = 0;

    if( nci==0.0 && mci==1.0 ) 
    {
        for( i=1; i<=XGRID; i++ )
            for( j=1; j<=YGRID; j++ )
                if( !boundary[i][j] )
                {
                    if( area[i][j]*cellarea>=chanThreshold && !flood[i][j] ) {
                        chan[i][j] = CHANNEL;
                        numChannelCells++;
                    }
                    else chan[i][j] = HILLSLOPE;
                } 
    }
    else
    {
        for( i=1; i<=XGRID; i++ )
            for( j=1; j<=YGRID; j++ )
                if( !boundary[i][j] )
                {
                    S = slope[i][j];
                    if( S>0 )
                    {
                        fn=pow(area[i][j]*cellarea,mci)*pow(S,nci);
                        if( fn > chanThreshold ) {
                            chan[i][j] = CHANNEL;
                            numChannelCells++;
                        }
                        else chan[i][j] = HILLSLOPE;
                    }
                    else chan[i][j] = HILLSLOPE;
                }
    }

} /**** end of FindChannels ****/ 



/************************************************************************
 *    
 *  CompareChannelAreas:    This function is used by the qsort() function
 *                          to determine which of two model cells has
 *                          the greater drainage area. It is identical
 *                          to CompareAreas, below, except that non-
 *                          channel cells are treated as having a 
 *                          large area, so that they end up at the end
 *                          of the sorted list and are ignored (since
 *                          the fluvial algorithm only counts up to
 *                          numChannelCells and ignores the rest).
 *                          For more information on the form of a
 *                          comparison function used by the qsort()
 *                          routine, consult the qsort manual page.
 *
 *      Calls: (none)
 *      Called by: qsort() library routine (in catchment-scale mode only)
 *
 ************************************************************************/

int CompareChannelAreas( a, b )
struct CellCoord *a, *b;
{
    int ax = a->x, ay = a->y, bx = b-> x, by = b->y;

    if( !chan[ax][ay] )
    {
        if( !chan[bx][by] ) return( 0 );
        else return( 1 );
    }
    else if( !chan[bx][by] ) return( -1 );

    if( area[ax][ay] > area[bx][by] ) return( 1 );
    else if( area[ax][ay] < area[bx][by] ) return( -1 );
    else return( 0 );
}


/************************************************************************
 *    
 *  CompareAreas:   This is a comparison function used by the library
 *                  routine qsort() to compare the values of two elements
 *                  in the list to be sorted. In this case, the value
 *                  compared is drainage area.
 *                      For more information on the form of a
 *                  comparison function used by the qsort()
 *                  routine, consult the qsort manual page.
 *
 *      Calls: (none)
 *      Called by: qsort() library routine
 *
 ************************************************************************/

int CompareAreas( a, b )
struct CellCoord *a, *b;
{
    if( area[a->x][a->y] > area[b->x][b->y] ) return( 1 );
    else if( area[a->x][a->y] < area[b->x][b->y] ) return( -1 );
    else return( 0 );
}


void shunt(int k,int m,struct CellCoord temp,struct CellCoord array[]){
  /*Heap building and re-ordering function for heapsort 
    Mike Bithell 03/10/00
  */
  int xjm1, yjm1, xj, yj;
  int i=k;
  int j=k+k;
  int xtmp = temp.x;
  int ytmp = temp.y;

  while (j <= m){
    xjm1=array[j-1].x;
    yjm1=array[j-1].y;
    if (j < m){
      xj=array[j].x;
      yj=array[j].y;
      if( area[xjm1][yjm1] < area[xj][yj] ) j=j+1;
      /*if(array[j-1] < array[j])j=j+1;*/
    }
    if( area[xtmp][ytmp] < area[xjm1][yjm1] ) {
      /*if (temp < array[j-1]){*/
      array[i-1]=array[j-1];
      i=j;
      j=j+j;
    }
    else{
      j=m+1;
    }
  }
  array[i-1]=temp;
}


void heapsort(int n,struct CellCoord array[]){
  /* A function that sorts the array, of dimension n, into increasing
     order using a heap sort (cribbed and adapted from numerical recipes)
     calls shunt to reorder the array elements. First set of calls orders
     the array into a heap and second set of calls does the sort.
     Performance should be O(NlogN) even for initially sorted arrays, but
     a little slower than quicksort.
     Mike Bithell 03/10/00
  */
  int k;
  struct CellCoord temp;

  if (n  < 2 ) return;

  printf("HEAPSORT HERE\n");
  
  for (k=n/2;k>=1;k--){
    temp=array[k-1];
    shunt(k,n,temp,array);
  }
  for(k=n;k>1;k--) {
    temp=array[k-1];
    array[k-1]=array[0];
    shunt(1,k-1,temp,array);
  }
  array[0]=temp;
}


/************************************************************************
 *    
 *  ErodeBedrock: This function computes the depth of stream erosion into
 *                bedrock over the duration of time step delt, as a 
 *                function of discharge, slope, and lithology. Negative
 *                erosion (ie, deposition due to negative slope term) is
 *                not allowed. The BEDROCK channel status is assigned to
 *                the current cell whenever this function is called.
 *                Version 5.10 (summer 1996) adds a "kbeff" term, which
 *
 *      Calls: (math library function: pow)
 *      Called by: Fluvial, UpdateElevations  
 *
 ************************************************************************/

double ErodeBedrock( x, y, h1, h2, delt )
int x, y;
double h1, h2, delt;
{
    double brterm, kq, S, qq;
    int lith = rockType[x][y];

    /*   qq = precip_flow*area[x][y];*/
    qq = disch[x][y];
    S = (h1-h2)/delx[x][y];
    if( S<0 ) return( 0 );
    brterm = layer[lith].K*(pow(qq,mb)*pow(S,nb)-layer[lith].tc)*delt;
    /*if( x==13 && y==4 ) printf("(%d,%d) A: %d Q: %f\tS:%f\tKb: %e\tE: %f\tDt: %.1f\tDh: %f\nchan: %d\tflood: %d\tdelx: %f\n",x,y,area[x][y],qq,S,layer[lith].K,brterm/delt,delt,brterm,chan[x][y],flood[x][y],delx[x][y]);*/
    chan[x][y] = BEDROCK;
    if( brterm<0 ) brterm=0;
    return( brterm );
}


/************************************************************************
 *                                                                      
 *  GenericTransport:  This a "generic" function for computing sediment
 *                     transport capacity, according to the equation
 *                     Qs = kf Q^m S^n. It does not include a threshold
 *                     term, though it may be modified to do so in the
 *                     future.
 *
 *      Calls: (math lib function: pow) 
 *      Called By: EstimateDtMax (via function pointer (*Transport)) 
 *                                                                      
/************************************************************************/

double GenericTransport( x, y, nx, ny )
int x, y, nx, ny;
{
    double S, Qw;

    S = (elev[x][y]-elev[nx][ny])/delx[x][y];
    /*    Qw = area[x][y]*precip_flow;*/
    Qw = disch[x][y];
    if( S <= 0.0 ) {
        return( 0.0 );
    }
    else return( kf*pow(Qw,mf)*pow(S,nf) );
}


/************************************************************************
 *                                                                      
 *  ThresholdTransport: Like GenericTransport() except with a threshold
 *                      term. Exponents may be chosen so that this
 *                      function becomes like the Meyer-Peter Mueller
 *                      and similar transport formulas, with
 *                      qs ~ (tau - taucrit)^p. The equation solved
 *                      includes channel width as the square root of
 *                      discharge (with the constant lumped into kf and
 *                      gkt),
 *
 *                      Qs = k W ( T - Tc )^p
 *
 *                      Note: because width dependence is explicitly
 *                      included in the leading term, the exponent mf
 *                      will not generally be the same as the exponent
 *                      used in the GenericTransport function.
 *
 *      Calls: (math lib function: pow) 
 *      Called By: EstimateDtMax (via function pointer (*Transport)) 
 *                                                                      
/************************************************************************/

double ThresholdTransport( x, y, nx, ny )
int x, y, nx, ny;
{
    double S, Qw, tauex;

    S = (elev[x][y]-elev[nx][ny])/delx[x][y];
    Qw = disch[x][y];
    S = S * (S>0);
    tau[x][y] = gkt*pow(Qw,mf)*pow(S,nf);
    /*Xtauex = gkt*pow(Qw,mf)*pow(S,nf) - gtaucrit;*/
    tauex = tau[x][y] - gtaucrit;
    tauex = (tauex>0) ? tauex : 0.0;
    return( kf*sqrt(Qw)*pow( tauex, gpf ) );
}

/************************************************************************
 *                                                                      
 *  MPMTransport
 *
 *  Computes transport rate using a version of the Meyer-Peter Mueller
 *  formula, assuming channel width proportional to the square root of
 *  mean annual discharge, which is the product of mean annual flow per
 *  unit area (gp_mean_ann) and area. All terms are assumed to be in
 *  SI units (unlike other transport functions in this code, which use
 *  time units of years!); the kf factor converts from m3/s to m/yr.cell.
 *
 *      Calls: (math lib function: pow) 
 *      Called By: EstimateDtMax (via function pointer (*Transport)) 
 *                                                                      
/************************************************************************/

double MPMTransport( x, y, nx, ny )
int x, y, nx, ny;
{
    double S, Qw, tauex;

    S = (elev[x][y]-elev[nx][ny])/delx[x][y];
    if( S<0.0 ) S=0.0;
    tau[x][y] = gkt*pow(area[x][y]*cellarea,0.3)*pow(S,0.7);
    tauex = tau[x][y] - gtaucrit;
    if( tauex<0.0 ) tauex = 0.0;
    return( kf*sqrt((double)area[x][y])*pow( tauex, 1.5 ) );
}


/************************************************************************
 *                                                                      
 *  StreamPwrDh:  Specialized fluvial sediment erosion/deposition 
 *                that finds a local analytical solution for the case in
 *                which Qs goes as Q x S (ie, linear). The accuracy of
 *                the local analytical solution remains to be fully
 *                tested.
 *
 *      Calls: (none) 
 *      Called By: Fluvial 
 *                                                                      
/************************************************************************/

double StreamPwrDh( x, y, nx, ny, l, delt )
int x, y, nx, ny;
double l, delt;
{
    double S, Qw;

    S = (elev[x][y]-elev[nx][ny])/delx[x][y];
    /*Qw = area[x][y]*precip_flow;*/
    Qw = disch[x][y];
    if( S <= 0.0 ) S = 0; 
    else return( (l-kf*Qw*S*delt)/(1.0+kf*Qw*(delt/delx[x][y])) );
}


/************************************************************************
 *                                                                      
 *  Bagnold:  Specialized function for sediment transport capacity that
 *            calculates capcity using the modified Bagnold bedload
 *            transport formula of Bridge & Dominic (1984). Shear stress
 *            is assumed to go as Q^(1/3) S^(2/3). Variable mf is 
 *            assumed to contain the critical shear stress.
 *
 *      Calls: (math library routines: pow, sqrt) 
 *      Called By: EstimateDtMax (via function pointer (*Transport) ) 
 *                                                                      
/************************************************************************/

double Bagnold( x, y, nx, ny )
int x, y, nx, ny;
{   
    double S, Qw, Qs, tauc=mf;

    S = (elev[x][y]-elev[nx][ny])/delx[x][y];
    if( S <= 0.0 ) return( 0.0 );
    /*Qw = area[x][y]*precip_flow;*/ 
    Qw = disch[x][y];
    tau[x][y] = pow(Qw,0.33)*pow(S,0.67);
    if( tau[x][y] > tauc )
        Qs = kf*sqrt(Qw)*(tau[x][y]-tauc)*(sqrt(tau[x][y])-sqrt(tauc));
    else Qs = 0.0;
    return( Qs );

}



double EstimateDTMax( qs )
double qs[XGRID+2][YGRID+2];
{
    int a,x,y,      /* Counters */
        nx,ny;      /* Neighbor cell coords */
    double ddt,
        locddt,
        S,
        kdt,    /* Should be an enterable param */
        dht[XGRID+2][YGRID+2],
        qsup[XGRID+2][YGRID+2];
    int temp;

    /* Initialize dht and qsup to zero */
    for( x=0; x<=XGRIDP1; x++ )
        for( y=0; y<=YGRIDP1; y++ )
        {
            dht[x][y] = 0.0;
            qsup[x][y] = 0.0;
        }

    /* Calculate dh/dt and Qs at start of current time step */
    for( a=0; a<numChannelCells; a++ )
    {
        x = cell[a].x;
        y = cell[a].y;
        nx = nbrx[x][y];
        ny = nbry[x][y];

        if( flood[x][y] )
            qs[x][y] = 0;
        else
            qs[x][y] = Transport( x, y, nx, ny );
        dht[x][y] = qsup[x][y] - qs[x][y];

        /* Correction for bedrock channels */
        if( dht[x][y]<0.0 && chansed[x][y]<=chanDepth[x][y] )
        {
            dht[x][y] = -ErodeBedrock( x,y,elev[x][y],elev[nx][ny],1.0 );
            /*qs[x][y] = qsup[x][y] - dht[x][y];*/
        }

        /* Send sediment flux to downstream node */
        qsup[nx][ny] += qs[x][y];
    }

    /* For hillslope cells, qs should be zero */
    for( a=numChannelCells; a<GRIDSIZE; a++ )
        qs[cell[a].x][cell[a].y] = 0.0;

    /* To improve numerical stability, take a weighted average of 
     * upstream + downstream (new in Version 5.10) */
    for( a=numChannelCells-1; a>=0; a-- )
    {
        x = cell[a].x;
        y = cell[a].y;
        nx = nbrx[x][y];
        ny = nbry[x][y];
        /* Clause added 3/98: don't smooth if one of the 2 cells is flooded */
        if( !flood[x][y] && !flood[nx][ny] )
        {
           dht[nx][ny] += downstrmWt*dht[x][y];
           qs[nx][ny] += downstrmWt*qs[x][y];
           dht[x][y] *= upstrmWt;
           qs[x][y] *= upstrmWt;
        }
    }

    /* Set dh/dt at east and west boundaries */
    for( y=0; y<=YGRIDP1; y++ )
    {
        dht[0][y] = 0.0;
        dht[XGRIDP1][y] = dht[XGRID][y];
    }
    
    /* For each cell, find DT that will bring slope to zero,*/
    /* or some fraction of the way toward zero. Take the    */ 
    /* smallest of these as minimum time step size.     */
    ddt = Dt;
    for( a=0; a<numChannelCells; a++ )
    {
        x = cell[a].x;
        y = cell[a].y;
        nx = nbrx[x][y];
        ny = nbry[x][y];
        

        if( dht[nx][ny] > dht[x][y] && !flood[x][y] )
        {
            locddt = 0.8*(elev[nx][ny]-elev[x][y])/(dht[x][y]-dht[nx][ny]);
            /*Xif( locddt < 0.0 ) locddt = (float)Dt;*/
            if( locddt < ddt ) {
                ddt = locddt;
            }
        }
    }

    return( ddt );

}



void UpdateElevations( delt, qs )
double delt, qs[XGRID+2][YGRID+2];
{
    int a,x,y,nx,ny;
    double dh,
        cap,
        depth,
        bedsed;  /* Depth of sediment on the bed (as opposed to the cell) */

    /* Sweep grid and reset loads for new computation */
    for( x=0; x<=XGRIDP1; x++ )
        for( y=0; y<=YGRIDP1; y++ )
            load[x][y] = 0;

    /* Added 1/95: assign incoming sed flux for discharge inputs (if any) */
    /* The 0.5 term comes from assuming that half of the mass denuded from */
    /* the upstream basin is removed via wash and dissolved load */
    for( a=0; a<nPtSources; a++ )
        load[XGRID][inputLoc[a]] = 0.5*upliftRate*(delt/Dt)*Ainput[a];

    for( a=0; a<numChannelCells; a++ )
    {
        x = cell[a].x;
        y = cell[a].y;
        nx = nbrx[x][y];
        ny = nbry[x][y];
        chan[x][y] = ALLUVIAL;

        dh = load[x][y] - qs[x][y]*delt;
        cap = dh;
        bedsed = (chansed[x][y]>chanDepth[x][y]) ?
            chansed[x][y] - chanDepth[x][y] : 0.0;
        if( -dh > bedsed )
        {
            dh = -bedsed
                - ErodeBedrock( x,y,elev[x][y],elev[nx][ny],delt);
            chan[x][y] = BEDROCK;
            if( dh < cap ) 
            {
                dh = cap;
                chan[x][y] = CAPLIM;
            }
        }

/*      if( a>10000 )
        {
            printf("%d %d = %.2f at dh=%f (Flood: %d) (chan: %d)\n",x,y,elev[x][y],dh,flood[x][y], chan[x][y]);
            printf("load: %f  qs: %f  delt: %f   a: %.2f\n",load[x][y],qs[x][y],delt,area[x][y]);
            printf("-> %d %d, h=%f\n",nx,ny,elev[nx][ny]);
            getchar();
        }
*/
        /*
         *  ADJUST ROCK AND SOIL THICKNESSES, SEND NEW 
         *  SEDIMENT LOAD
         *  DOWNSTREAM, AND INCREMENT DENUDATION COUNTER
         */
        elev[x][y] += dh;

        if( elev[x][y]<-100.0 )
        {
            printf("ALERT %d %d hits %f at dh=%f\n",x,y,elev[x][y],dh);
            printf("load: %f  qs: %f  delt: %f   a: %d\n",load[x][y],qs[x][y],delt,area[x][y]);
            printf("-> %d %d, h=%f\n",nx,ny,elev[nx][ny]);
            printf("chan: %d    flood:  %d    lithology: %d, chansed=%f\n",chan[x][y],flood[x][y],rockType[x][y],chansed[x][y]);
            printf( "Bedrock erosion term: %f\n",dh+chansed[x][y]);
            problem( "Negative elevation in UpdateElevation()\n");
        }
        delz[x][y][0] -= (float)dh;
        denud -= dh;
        load[nx][ny] += load[x][y] - dh;
    
        /*
         * ADD DEPOSITED MATERIAL TO SOIL, OR REMOVE MATERIAL
         * FROM SOIL AND UNDERLYING ROCK LAYERS, AS APPROP.
         */
        if( chansed[x][y] >= -dh )
            chansed[x][y] += dh;
        else
        {
            depth = -dh - chansed[x][y];
            chansed[x][y] = 0.0;       
            if( rock[x][y] >= depth )
                rock[x][y] -= depth;

            /*
             *  IN THIS CASE, BOTH SOIL AND UPPER ROCK
             *  LAYER ARE REMOVED, AND THE UNDERLYING 
             *  LAYER NOW BECOMES THE UPPER LAYER.
             *  (MODEL ASSUMES THAT EROSION WILL NEVER
             *  PENETRATE MORE THAN ONE ROCK LAYER IN
             *  A TIMESTEP).
             */
            else
            {
                depth -= rock[x][y];
                rockType[x][y]++;
                rock[x][y] = layer[rockType[x][y]].T - depth;
            }
        }
                        
    }

    /* In Catchmode mode, there may be discontinuous channel segments.
     * To make sure mass is conserved, we deposit any sediment recorded
     * in the LOAD variable.
     */
    if( mode==CatchmentMode )
        for( a=numChannelCells; a<GRIDSIZE; a++ )
        {
            x = cell[a].x;
            y = cell[a].y;
            if( !boundary[x][y] )
            {
                elev[x][y] += load[x][y];
                chansed[x][y] += load[x][y];
            }
        }

    /* Add up all the sed that's been delivered to the west boundary */
    for( a=1; a<=YGRID; a++ )
    {
        sed_yield += load[0][a];
        sed_yield += load[XGRIDP1][a];
        load[0][a] = 0;
        load[XGRIDP1][a] = 0;
    }

}



double Solve( delt, dtm, qs, level )
double delt, dtm, qs[XGRID+2][YGRID+2];
int level;
{
    int i;

    /*if( level > 20 ) 
        problem("Infinite recursion in SOLVE\nCheck POND setting\n");

    if( dtMinApplied > delt ) dtMinApplied = delt;*/

    if( delt > dtm && delt > dt_min )
    {
        dtm = Solve( 0.5*delt, dtm, qs, level+1 );
        dtm = Solve( 0.5*delt, dtm, qs, level+1 );
    }
    else
    {
        UpdateElevations( delt, qs );
        if( level>0 ) dtm = EstimateDTMax( qs );
    }
    return( dtm );

}



Fluvial()
{
    int a,          /* Counter for nodes on the list */
        x,      /* x coord of current cell */
        y,      /* y coord of current cell */
        i,      /* counter for sub-steps in bedrock-only mode */
        j,k,
        nSubSteps, /* No. of sub-time steps (used only in STREAMPWR mode */
        nx,     /* x-coord of neighbor node     */
        ny;     /* y-coord of neighbor node     */
    double qs[XGRID+2][YGRID+2],
        dhdt[XGRID+2][YGRID+2],  /* Bedrock erosion rates */
        cap,
        depth,
        beta,
        dh,
        rt,    /* Remaining time within 'global' time-step */
        locdt, /* Local (at a cell) maximum time-step size */
        dtmax;

    /*printf("Fluvial()...");*/
    
    if( kf==0 )
    {
       rt = Dt;
       do
       {
          /* For each node, compute erosion rate and time to zero slope */
          /* We can do this in one loop if we process in downstream-to- */
          /* upstream order. */
          dtmax = rt;
          for( a=numChannelCells-1; a>=0; a-- )
          {
             /* Compute the bedrock erosion rate */
             x = cell[a].x;
             y = cell[a].y;
             nx = nbrx[x][y];
             ny = nbry[x][y];
             dhdt[x][y] = -ErodeBedrock( x,y,elev[x][y],elev[nx][ny],1.0 );
             
             /* If downstream node is eroding less rapidly, figure out how */
             /* long it would take to reach a zero slope, and take a       */
             /* fraction of that time as local maximum.                    */
             if( dhdt[nx][ny] > dhdt[x][y] )
             {
                locdt = 0.2 * ( ( elev[x][y] - elev[nx][ny] )
                                / ( dhdt[nx][ny] - dhdt[x][y] ) );
                if( locdt < dtmax ) dtmax = locdt;
             }
          }
          if( dtmax<dt_min ) dtmax = dt_min;
          /*printf( "dt=%f\n",dtmax);*/
          
          /* Now extrapolate the erosion rates to compute new elevations */
          for( a=numChannelCells-1; a>=0; a-- )
          {
             x = cell[a].x;
             y = cell[a].y;
             nx = nbrx[x][y];
             ny = nbry[x][y];
             dh = dhdt[x][y] * dtmax;
             
             /* Since we're in bedrock-only mode, remove any sed */
             if( chansed[x][y]>0 ) {
                elev[x][y] -= chansed[x][y];
                chansed[x][y] = 0;
             }
             
             /* Make sure slope doesn't drop below minimum */
             /* (this is the only reason to work upstream in this loop) */
             if( optMinSlope && elev[x][y]<elev[nx][ny]+delx[x][y]*minS )
             {
                dh = (elev[nx][ny]+delx[x][y]*minS)-elev[x][y];
                if( dh > 0 ) dh=0;
             }
             
             /* Update elevation, rock thickness and rock type */
             elev[x][y] += dh;
             rock[x][y] += dh;
             if( rock[x][y] < 0 ) {
                rockType[x][y]++;
                rock[x][y] = layer[rockType[x][y]].T;
             }
          }
          rt -= dtmax;
          
       } while( rt>0.0 );

    } /* if kf=0 (bedrock-mode) */
    
    else if( optTransFn==STREAMPWR )
    {
        nSubSteps = (int)((float)Dt/dt_min);
        for( i=1; i<=nSubSteps; i++ )
        {
            for( j=0; j<=XGRIDP1; j++ ) for( k=1; k<=YGRID; k++ )
                load[j][k] = 0;
            for( a=0; a<numChannelCells; a++ )
            {
                x = cell[a].x;
                y = cell[a].y;
                nx = nbrx[x][y];
                ny = nbry[x][y];
                dh = StreamPwrDh(x,y,nx,ny,load[x][y],dt_min);
                if( -dh > chansed[x][y] )
                {
                    cap = dh;
                    dh = -chansed[x][y] - ErodeBedrock( x,y,
                          elev[x][y],elev[nx][ny],dt_min);
                    if( dh < cap ) 
                    {
                        dh = cap;
                        chan[x][y] = CAPLIM;
                    }
                }
                elev[x][y] += dh;
                delz[x][y][0] -= (float)dh;
                denud -= dh;
                load[nx][ny] += load[x][y] - dh;
    
                /*
                 * ADD DEPOSITED MATERIAL TO SOIL, OR REMOVE MATERIAL
                 * FROM SOIL AND UNDERLYING ROCK LAYERS, AS APPROP.
                 */
                if( chansed[x][y] >= -dh )
                    chansed[x][y] += dh;
                else
                {
                    depth = -dh - chansed[x][y];
                    chansed[x][y] = 0.0;       
                    if( rock[x][y] >= depth )
                        rock[x][y] -= depth;
                    else
                    {
                        depth -= rock[x][y];
                        rockType[x][y]++;
                        rock[x][y] = layer[rockType[x][y]].T - depth;
                    }
                }
                        
            }
        }
        
        /* Add up the sed flux along the left boundary */   
        sed_yield = 0;
        for( i=1; i<=YGRID; i++ )
            sed_yield += load[0][i];
    }
    else
    {
        dtMinApplied = (float)Dt;
        dtmax = EstimateDTMax( qs );
        Solve( (float)Dt, dtmax, qs, 0 );
    }

    /*printf("done\n");*/

} /* Fluvial */


/************************************************************************
 *    
 * CalcPrecip: This function computes sinusoidally varying or step-function
 *             varying precipitation. It is called only when the flag 
 *             optVarPrecip is set to something other than zero. Note
 *             that precipitation rate here refers to the intensity, 
 *             rather than the mean. Mean annual runoff is assumed to
 *             stay constant, so that if precip intensity goes up, the
 *             frequency or duration goes down in proportion. Since
 *             frequency/duration are reflected in the fluvial transport
 *             and erosion coefficients, these coefficients are adjusted
 *             accordingly. Function added for version 5.10 (was prev'ly
 *             in "special" version WE38-5).
 *                 Global parameters used are:
 *                   precip_mean:  mean precipitation rate (m/yr)
 *                   precip_var: maximum change in precipitation
 *                               during cycle (m/yr); max precip = 
 *                               precip_mean + precip_var; min =
 *                               precip_mean - precip_var
 *                   lambdap: period of oscillation (yrs)
 *                   kfmean: mean transport efficiency coefficient
 *                   kf: transport coefficient adjusted for intermittency
 *                   layer[].K and .Kmn: mean and adjusted bedrock
 *                                       erodibilities
 *
 *        Calls: (math library function: fmod)
 *        Called by: main
 *
 ************************************************************************/
CalcPrecip( tstep )
int tstep;
{
    float time = tstep*Dt;
    int i;

    /* If step fn, calculate whether it's "up" or "down", ... */
    if( optVarPrecip==STEPWAVE )
    {
      /* TODO: shouldn't the statement below use time not tstep? */
        if( fmod( tstep, lambdap )< 0.5*lambdap )
            precip_flow = precip_mean+precip_var;
        else
            precip_flow = precip_mean-precip_var;
    }
    else  /* ...otherwise, sinusoidally varying */
        precip_flow = precip_mean+precip_var*sin((TWOPI*time)/lambdap);

    /* Adjust intermittency factor(s) */
    kf = kfmean*(precip_mean/precip_flow);
    for( i=1; i<=n_layers; i++ )
        layer[i].K = layer[i].Kmn*(precip_mean/precip_flow);

}



void DiffuseExpl()
{
    double dhin[XGRID+2][YGRID+2],
        dhout[XGRIDP1][YGRIDP1],
        qsn,qss,qse,qsw,
        dhtot,
        scale,k;
    int i,j,jn,js;

    k = kd*Dt/cellarea;

    for( i=1; i<=XGRID; i++ )
        for( j=1; j<=YGRID; j++ )
            dhin[i][j]=0.0;
    for( i=1; i<=XGRID; i++ )
        for( j=1; j<=YGRID; j++ )
        {
            scale=1.0;
            jn=(j==YGRID) ? 1 : j+1;
            js=(j==1) ? YGRID : j-1;
            qsn= (elev[i][j]>elev[i][jn]) ? k*(elev[i][j]-elev[i][jn]) : 0.0;
            qss= (elev[i][j]>elev[i][js]) ? k*(elev[i][j]-elev[i][js]) : 0.0;
            qse= (elev[i][j]>elev[i+1][j]) ? k*(elev[i][j]-elev[i+1][j]) : 0.0;
            qsw= (elev[i][j]>elev[i-1][j]) ? k*(elev[i][j]-elev[i-1][j]) : 0.0;
            if((dhtot=qsn+qss+qse+qsw)>chansed[i][j]) scale=chansed[i][j]/dhtot;
            dhout[i][j]=scale*dhtot;
            dhin[i][jn]+=scale*qsn;
            dhin[i][js]+=scale*qss;
            dhin[i+1][j]+=scale*qse;
            dhin[i-1][j]+=scale*qsw;

        }

    /* Adjust elevations, sediment reservoirs and denudation rates */
    for( j=1; j<=YGRID; j++ )
    {
        for( i=1; i<=XGRID; i++ )
        {
            chansed[i][j] += dhin[i][j]-dhout[i][j];
            elev[i][j] += dhin[i][j]-dhout[i][j];
            delz[i][j][0] += (float)(dhout[i][j]-dhin[i][j]);
        }

        /* Accumulate influx at boundaries in sed_yield vbl */
        sed_yield += dhin[0][j] + dhin[YGRIDP1][j];
    }

}


/*
 *  DiffuseExplIrreg: Version of DiffuseExpl designed for irregular boundary.
 *      Doesn't compute a sediment flux for boundary cells, or from grid
 *      cells to boundary cells.
 */
void DiffuseExplIrreg()
{
    double dhin[XGRID+2][YGRID+2],
        dhout[XGRIDP1][YGRIDP1],
        qsn,qss,qse,qsw,
        dhtot,
        scale,k;
    int i,j,jn,js;

    k = kd*Dt/cellarea;

    for( i=1; i<=XGRID; i++ )
        for( j=1; j<=YGRID; j++ )
            dhin[i][j]=0.0;
    for( i=1; i<=XGRID; i++ )
        for( j=1; j<=YGRID; j++ )
            if( !boundary[i][j] )
            {
                qsn=qss=qse=qsw=0.0;
                scale=1.0;
                jn=(j==YGRID) ? 1 : j+1;
                js=(j==1) ? YGRID : j-1;
                if( boundary[i][jn]!=NoFluxBoundary )
                    qsn= (elev[i][j]>elev[i][jn]) ? k*(elev[i][j]-elev[i][jn]) : 0.0;
                if( boundary[i][js]!=NoFluxBoundary )
                    qss= (elev[i][j]>elev[i][js]) ? k*(elev[i][j]-elev[i][js]) : 0.0;
                if( boundary[i+1][j]!=NoFluxBoundary )
                    qse= (elev[i][j]>elev[i+1][j]) ? k*(elev[i][j]-elev[i+1][j]) : 0.0;
                if( boundary[i-1][j]!=NoFluxBoundary )
                    qsw= (elev[i][j]>elev[i-1][j]) ? k*(elev[i][j]-elev[i-1][j]) : 0.0;
                /*Xif((dhtot=qsn+qss+qse+qsw)>chansed[i][j]) scale=chansed[i][j]/dhtot;
                dhout[i][j]=scale*dhtot;
                dhin[i][jn]+=scale*qsn;
                dhin[i][js]+=scale*qss;
                dhin[i+1][j]+=scale*qse;
                dhin[i-1][j]+=scale*qsw;*/
                
                dhtot = qsn+qss+qse+qsw;
                if( optDiffuse==WEATHLIMDIFF )
                {
                    if( dhtot>chansed[i][j] ) scale=chansed[i][j]/dhtot;
                    dhout[i][j]=scale*dhtot;
                    dhin[i][jn]+=scale*qsn;
                    dhin[i][js]+=scale*qss;
                    dhin[i+1][j]+=scale*qse;
                    dhin[i-1][j]+=scale*qsw;
                }
                else
                {
                    dhout[i][j]=dhtot;
                    dhin[i][jn]+=qsn;
                    dhin[i][js]+=qss;
                    dhin[i+1][j]+=qse;
                    dhin[i-1][j]+=qsw;
                }
            }

    /* Adjust elevations, sediment reservoirs and denudation rates */
    for( j=1; j<=YGRID; j++ )
    {
        for( i=1; i<=XGRID; i++ )
            if( !boundary[i][j] )   
            {
                chansed[i][j] += dhin[i][j]-dhout[i][j];
                if( chansed[i][j]<0 )
                  {
                    rock[i][j]+=chansed[i][j];
                    chansed[i][j]=0;
                    if( rock[i][j]<0 ) 
                      {
                        rockType[i][j]++;
                        rock[i][j] += layer[rockType[i][j]].T;
                      }
                  }
                elev[i][j] += dhin[i][j]-dhout[i][j];
                delz[i][j][0] += (float)(dhout[i][j]-dhin[i][j]);
            }

        /* Accumulate influx at boundaries in sed_yield vbl */
        /*sed_yield += dhin[0][j] + dhin[YGRIDP1][j];*/
    }

}



void DiffuseAll()
{
    double dhin[XGRID+2][YGRID+2],
        dhout[XGRIDP1][YGRIDP1],
        qsn,qss,qse,qsw,
        dhtot,
        k;
    int i,j,jn,js;

    k = kd*Dt/(dx*dx);

    for( i=0; i<=XGRIDP1; i++ )
        for( j=0; j<=YGRIDP1; j++ )
            dhin[i][j]=0.0;
    for( i=1; i<=XGRID; i++ )
        for( j=1; j<=YGRID; j++ )
            if( !boundary[i][j] )
            {
                qsn=qss=qse=qsw=0.0;
                jn=(j==YGRID) ? 1 : j+1;
                js=(j==1) ? YGRID : j-1;
                if( boundary[i][jn]!=NoFluxBoundary )
                    qsn= (elev[i][j]>elev[i][jn]) ? k*(elev[i][j]-elev[i][jn]) : 0.0;
                if( boundary[i][js]!=NoFluxBoundary )
                    qss= (elev[i][j]>elev[i][js]) ? k*(elev[i][j]-elev[i][js]) : 0.0;
                if( boundary[i+1][j]!=NoFluxBoundary )
                    qse= (elev[i][j]>elev[i+1][j]) ? k*(elev[i][j]-elev[i+1][j]) : 0.0;
                if( boundary[i-1][j]!=NoFluxBoundary )
                    qsw= (elev[i][j]>elev[i-1][j]) ? k*(elev[i][j]-elev[i-1][j]) : 0.0;
                dhtot=qsn+qss+qse+qsw;
                dhout[i][j]=dhtot;
                dhin[i][jn]+=qsn;
                dhin[i][js]+=qss;
                dhin[i+1][j]+=qse;
                dhin[i-1][j]+=qsw;

            }

    /* Adjust elevations, sediment reservoirs and denudation rates */
    for( j=1; j<=YGRID; j++ )
    {
        for( i=1; i<=XGRID; i++ )
            if( !boundary[i][j] )   
            {
                chansed[i][j] += dhin[i][j]-dhout[i][j];
                if( chansed[i][j]<0 )
                  {
                    rock[i][j]+=chansed[i][j];
                    chansed[i][j]=0;
                    if( rock[i][j]<0 ) 
                      {
                        rockType[i][j]++;
                        rock[i][j] += layer[rockType[i][j]].T;
                      }
                  }
                elev[i][j] += dhin[i][j]-dhout[i][j];
                delz[i][j][0] += (float)(dhout[i][j]-dhin[i][j]);
            }

        /* Accumulate influx at boundaries in sed_yield vbl */
        sed_yield += dhin[0][j] + dhin[XGRIDP1][j];
    }

}


/*
 *  DiffuseADI()
 * 
 *  Solves hillslope diffusion equation using Alternating-Direction Implicit
 *  method. Not currently operational. (TODO)
 */
void DiffuseADI()
{
        double R,
                P,              /* Scaling factor for kdusivity coef.s */
                uu[GRIDSIZE+1],
                rr[GRIDSIZE+1],
                rhs[XGRIDP1][YGRIDP1],
                u[XGRID+2][YGRID+2],
                Rew[XGRIDP1][YGRIDP1],  /* Diffusivities=f(soil) */
                Rns[XGRIDP1][YGRIDP1];  /* Diffusivities=f(soil) */
        int r, c, k,            /* Counters (r=rows, c=columns) */
                M = YGRID,      /* M = no. of rows */
                N = XGRID;      /* N = no. of columns */

        /* TODO: Initialize R properly */
        
        /* WRAP NORTH AND SOUTH BOUNDARIES */
        wrapNSBoundaries();

        /* COMPUTE RHS FOR ODD TRAVERSE */
        for( r=1; r<=M; r++ ) {
                for( c=1; c<=N; c++ ) 
                        rhs[c][r] = elev[c][r-1] 
                                        + (1.0/R-2.0)*elev[c][r]
                                        + elev[c][r+1];
                rhs[1][r] += elev[0][r];              /* WEST */
                rhs[N][r] += elev[N+1][r];            /* EAST */
        }

        /* Translate RHS matrix to vector, for odd traverse */
        k=1;
        for( r=1; r<=M; r++ )
                for( c=1; c<=N; c++ ) {
                        rr[k] = rhs[c][r];
                        k++;
                }

        /* Solve for odd traverse */
        tridag( uca, ucb, ucc, rr, uu, GRIDSIZE );

        /* Translate odd solution vector back to matrix */
        k=1;
        for( r=1; r<=M; r++ )
                for( c=1; c<=N; c++ ) {
                        u[c][r] = uu[k];
                        k++;
                }

        /* Adjust any floating boundaries */
        /* here, top and bottom bds will float */
        for( c=1; c<=N; c++ ) {
                /*u[c][1] = u[c][2];     Top */
                u[c][0] = u[c][YGRID];  /* Top */
                /*u[c][M] = u[c][M-1];   Bottom */
                u[c][M+1] = u[c][1];    /* Bottom */
        }
        
        /*
         *  SET EAST AND WEST BOUNDARIES. WEST BOUNDARY IS NORMALLY FIXED;
         *  TODO: Update boundary handling here.
         */
        for( r=0; r<=M+1; r++ ) {
                u[0][r] = elev[0][r];
                if( east_boundary==NO_FLUX )
                        u[N+1][r] = u[N][r]+0.5*(u[N][r]-u[N-2][r]);
                else
                        u[N+1][r] = elev[N+1][r];
        }
        
        /* Compute RHS for even traverse */
        for( c=1; c<=N; c++ ) 
        { 
                for( r=1; r<=M; r++ )
                        rhs[c][r] = u[c-1][r]
                                        + (1.0/R-2.0)*u[c][r]
                                        + u[c+1][r];
                rhs[c][1] += elev[c][0];
                rhs[c][M] += elev[c][M+1];
        }

        /* Translate RHS matrix to vector form for even traverse */
        k=1;
        for( c=1; c<=N; c++ )
                for( r=1; r<=M; r++ ) {
                        rr[k] = rhs[c][r];
                        k++;
                }

        /* SOLVE FOR NEW VALUES */
        tridag( vca, vcb, vcc, rr, uu, GRIDSIZE );

        /* Translate solution vector back to matrix form */
        k=1;
        for( c=1; c<=N; c++ )
                for( r=1; r<=M; r++ ) {
                        u[c][r] = uu[k];
                        k++;
                }

        /* New elevations (west boundary remains unchanged) */
        for( c=1; c<=XGRID; c++ )
                for( r=1; r<=YGRID; r++ ) 
                {
                        chansed[c][r] += u[c][r] - elev[c][r];
                        denud += elev[c][r] - u[c][r];
                        delz[c][r][0] += elev[c][r] - u[c][r];
                        if( chansed[c][r]<0.0 )
                        {
                                rock[c][r] += chansed[c][r];
                                chansed[c][r] = 0.0;
                                while( rock[c][r] < 0.0 )
                                {
                                        rockType[c][r]++;
                                        rock[c][r] = layer[rockType[c][r]].T
                                                        + rock[c][r];
                                }
                        }
                }


} /*  */


/* tridag: tridiagonal matrix solver */

static void tridag(a,b,c,r,u,n)
double *a,*b,*c,*r,*u;
int n;
{
        int j;
        double bet,gam[GRIDSIZE+1];

        if(b[1]==0.0) problem("Error 1 in TRIDAG.");
        u[1]=r[1]/(bet=b[1]);
        for(j=2;j<=n;j++) {
                gam[j]=c[j-1]/bet;
                bet=b[j]-a[j]*gam[j];
                if(bet==0.0) problem("Error 2 in TRIDAG." );
                u[j]=(r[j]-a[j]*u[j-1])/bet;
        }
        for(j=(n-1);j>=1;j--)
                u[j] -= gam[j+1]*u[j+1];

} /* tridag */




int CompareElevations( a, b )
struct CellCoord *a, *b;
{
    int ax = a->x, ay = a->y, bx = b-> x, by = b->y;

    if( elev[ax][ay] > elev[bx][by] ) return( 1 );
    else if( elev[ax][ay] < elev[bx][by] ) return( -1 );
    else return( 0 );
}



int checkSlopes( oversteep, type, hlist )
int oversteep[XGRID+2][YGRID+2], type;
struct CellCoord hlist[GRIDSIZE];
{
    double /*X hh,*/ 
      hcr,         /* Critical height of cell for failure */
      hcrs,        /* Critical height, non-diagonal slope orientation */
      hcrd,        /* Critical height, diagonal slope orientation */
      heff;           /* Effective elev, either of surface
                       or of soil-bedrock contact */
    int i, j, nsteep;

    if( type==SOIL ) hcrs=soil_slope*dx;
    hcrd = hcrs*ROOT2;

    /* Check and flag any oversteepened slopes */
    nsteep = 0;
    for( i=1; i<=XGRID; i++ )
        for( j=1; j<=YGRID; j++ )
        {
            if( type==ROCK ) {
                hcrs=layer[rockType[i][j]].Scr*dx;
                hcrd=hcrs*ROOT2;
            }
            heff=(type==ROCK) ? elev[i][j]-chansed[i][j] : elev[i][j];
            oversteep[i][j] = FALSE;
            if( (type !=SOIL || chansed[i][j]>0.0) && !flood[i][j] )
            {
	      /*X hh = 0.0;*/
                hcr = ( i!=nbrx[i][j] && j!=nbry[i][j] ) ? hcrd : hcrs;
                if( (heff - elev[nbrx[i][j]][nbry[i][j]]) > hcr )
                {
                    oversteep[i][j] = TRUE;
                    hlist[nsteep].x = i;
                    hlist[nsteep].y = j;
                    nsteep++;
                    slp_fail[i][j] = type;
                }
            }
        }

    /* Sort by elev */
    if( nsteep > 0 ) {
        qsort( hlist, nsteep, sizeof( hlist[0] ), CompareElevations );
    }
    return( nsteep );
}



int checkSlopesAreaDep( oversteep, type, hlist )
int oversteep[XGRID+2][YGRID+2], type;
struct CellCoord hlist[GRIDSIZE];
{
    double slope, acr;
    int i, j, nsteep;

    nsteep = 0;
    for( i=1; i<=XGRID; i++ )
        for( j=1; j<=YGRID; j++ )
        {
            oversteep[i][j] = FALSE;
            slope = (elev[i][j]-elev[nbrx[i][j]][nbry[i][j]])/delx[i][j];
            if( slope>angleOfRepose && chansed[i][j]>0 && !flood[i][j] )
            {
                acr = kls*sin(atan(slope))*(1.0-slope/soil_slope);
                if( area[i][j]*cellarea>acr )
                {
                    oversteep[i][j] = TRUE;
                    hlist[nsteep].x = i;
                    hlist[nsteep].y = j;
                    nsteep++;
                    slp_fail[i][j] = type;
                }
            }
        }

    /* Sort by elev */
    if( nsteep > 0 ) {
        qsort( hlist, nsteep, sizeof( hlist[0] ), CompareElevations );
    }
    return( nsteep );

}


void Collapse( nsteep, oversteep, type, hlist )
int nsteep, oversteep[XGRID+2][YGRID+2], type;
struct CellCoord hlist[GRIDSIZE];
{
    double h1, h2, hcr, dh1, dh2,
        hcrs,     /* Critical height, straight (non-diag) orientation */ 
        hcrd,     /* Critical height, diagonal orientation */
        hr,
        heff;           /* "Effective" height, either surface
                        elev or sed-rock interface elev */
    int x, y, done=FALSE, n,
                nextx, nexty, stable, bound, iter=0;


    if( type==SOIL ) hcrs=angleOfRepose*dx;
    hcrd=hcrs*ROOT2;
    for( n=0; n<nsteep; n++ )
    {
        bound = FALSE;
        x = hlist[n].x;
        y = hlist[n].y;
        if( type==ROCK ) {
            hcrs=layer[rockType[x][y]].Scr*dx;
            hcrd=hcrs*ROOT2;
        }
        nextx = nbrx[x][y];
        nexty = nbry[x][y];
        hcr = (x!=nextx && y!=nexty ) ? hcrd : hcrs;
        heff=(type==ROCK) ? elev[x][y]-chansed[x][y] : elev[x][y];
        dh1 = ((heff-elev[nextx][nexty])-hcr);
        if( !boundary[nextx][nexty] ) dh1 *= 0.5;  /* changed from isperim for 5.10, 8/96 */
        else bound = TRUE;
        if( type==SOIL && dh1 > chansed[x][y] ) dh1=chansed[x][y];
        /*if( dh1>0 ) printf("%d, %d COLLAPSE: %d\n",x,y,type);*/
        if( dh1 < VERYSMALL ) dh1=VERYSMALL; /* Prevents endless loops */
        elev[x][y] -= dh1;
        delz[x][y][0] += dh1;
        if( type==ROCK ) rock[x][y] -= dh1;
        if( rock[x][y] <= 0.0 )
        {
            rockType[x][y]++;
            rock[x][y] += layer[rockType[x][y]].T;
        }
        if( type==SOIL ) chansed[x][y] -=dh1;
        stable = FALSE;
        while( !stable && !bound )
        {
            x = nextx;
            y = nexty;
            nextx = nbrx[x][y];
            nexty = nbry[x][y];
            if( slope[x][y] <= 0.0 )
            {
                elev[x][y] += dh1;
                FindNextCell( x, y );
                elev[x][y] -= dh1;
                nextx = nbrx[x][y];
                nexty = nbry[x][y];
            }
            hr = (x!=nextx && y!=nexty ) ? hcrd : hcrs;
            dh2=(elev[x][y]+dh1-elev[nextx][nexty]-hr);
            dh2 *= 0.5;
            if( boundary[nextx][nexty] ) {   /* changed from isperim, 8/96 */
                dh2 *= 2;
                bound=TRUE;
            }
            else if( oversteep[nextx][nexty] )
                dh2 *= 2;

            if( dh2 <= 0.0 || flood[x][y] ) 
            {
                elev[x][y] += dh1;
                delz[x][y][0] -= dh1;
                chansed[x][y] += dh1;
                stable=TRUE;
            }
            else if( dh1 > dh2 )
            {
                elev[x][y] += dh1 - dh2;
                delz[x][y][0] += dh2 - dh1;
                chansed[x][y] += dh1 - dh2;
                dh1 = dh2;
            }
        }
    }
    iter++;

}


/*
**  CollapseSimple:  Fast slope failure algorithm that doesn't track sediment
**                   mass. Used in "bedrock-only" mode where all sediment
**                   evaporates anyway.
*/ 
void CollapseSimple()
{
    int a, x, y, nx, ny;
    float slp, dh, critdrop;

    /* From lowest to highest... */
    for( a=numChannelCells-1; a>=0; a-- )
    {
        x = cell[a].x;
        y = cell[a].y;
        nx = nbrx[x][y];
        ny = nbry[x][y];
        slp = (elev[x][y]-elev[nx][ny])/delx[x][y];
        if( slp > layer[rockType[x][y]].Scr )
        {
            critdrop = layer[rockType[x][y]].Scr*delx[x][y];
            dh = (elev[x][y]-elev[nx][ny])-critdrop;  /* excess height */
            if( dh<0 ) printf("WARNING: problem in CollapseSimple\n");
            elev[x][y] -= dh;
            rock[x][y] -= dh;
        }
    }

}


void SlopeCollapse()
{
    int i, j, oversteep[XGRID+2][YGRID+2];
    struct CellCoord hlist[GRIDSIZE];
    int ssteep, rsteep;
    int didCollapse=FALSE;

    /*printf("SlopeCollapse()...");*/

    /* Reset the slp_fail array */
    for( i=1; i<=XGRID; i++ )
        for( j=1; j<=YGRID; j++ )
            slp_fail[i][j]=FALSE;

    if( kf<=0 ) {  /* In this case, everything is bedrock---it's simple */
        CollapseSimple();
        return;
    }

    do {
        if( optAreaDepLS )
            ssteep = checkSlopesAreaDep( oversteep, SOIL, hlist );
        else 
            ssteep = checkSlopes( oversteep, SOIL, hlist );
        if( ssteep ) {
           Collapse( ssteep, oversteep, SOIL, hlist );
           didCollapse = TRUE;
           for( i=1; i<=XGRID; i++ )
                for( j=1; j<=YGRID; j++ )
                    FindNextCell(i,j);
        }
        rsteep = checkSlopes( oversteep, ROCK, hlist );
        if( rsteep ) {
           Collapse( rsteep, oversteep, ROCK, hlist );
           didCollapse = TRUE;
           for( i=1; i<=XGRID; i++ )
               for( j=1; j<=YGRID; j++ )
                   FindNextCell(i,j);
        }
    } while( rsteep>0 || ssteep>0 );

    /* Bug fix 3/98: if we did have an avalanche, any ponding information */
    /* will have been erased by the call(s) to FindNextCell(), so redo */
    if( pond )
        for( i=1; i<=XGRID; i++ )
            for( j=1; j<=YGRID; j++ )
                if( !boundary[i][j] && flood[i][j]==BASIN )
                    DoPonds(i,j);

    /*printf("done\n");*/
    
}



float getGrabenFactor( x )
int x;
{
        int block;

        block = (int)((x+grabenLength-1)/(float)grabenLength);
        if( block > ngrabens ) block=ngrabens;
        return( block/(float)ngrabens );
}


Uplift( uptype, t )
int uptype,t;
{
        if( uptype==PLATEAU )
                plateauUplift();
        else if( uptype==BLOCK )
                blockUplift( t );
        else if( uptype==TILTBLOCK )
                tiltBlockUplift();
        else if( uptype==DYNAMIC )
                dynamicUplift();
        else if( uptype==ERFSUB )
                erfSubsidence( t );
        else if( uptype==ERRFUNC )
                erfUplift();
        else if( uptype==HALFDOME )
                halfdomeUplift();
        else if( uptype==FOLDS )
            foldUplift( t, 0, 0 );
        else if( uptype==STRIKESLIP )
            StrikeSlip();
        else            
                problem("Specified uplift type not recognized.");       

}


plateauUplift()
{
        int i, j;
        
        for(i=1; i<=XGRID; i++ )
                for( j=0;j<=YGRIDP1; j++ )
                        if( !boundary[i][j] )
                            base[i][j] += upliftRate;

}


/*      Identical to plateau uplift, except that only west boundary     */
/*              remains fixed.                                          */
/*      Modified 6/9/94 to allow for diachronous uplift, ie, migration  */
/*          of the boundary fault through time. Fault remains fixed when*/
/*          vf is zero.                                                 */

blockUplift( t )
int t;
{
    int i, j, xf;

    xf = (int)(((float)fault*dx - vf*t*Dt)/dx);
    if( xf<1 ) xf=1;

    for(i=xf; i<=XGRIDP1; i++ )
        for( j=0; j<=YGRIDP1; j++ )
            base[i][j] += upliftRate;

}


tiltBlockUplift()
{
        int i, j;
        double d, dh;

        for(i=fault; i<=XGRIDP1; i++ )
                for( j=0; j<=YGRIDP1; j++ )
                {
                        d = fault_to_pivot - (i-fault)*dx;
                        dh = upliftRate*d/fault_to_pivot;
                        base[i][j] += dh;
                }

} /* tiltBlockUplift */



dynamicUplift()
{
        int i, j;
        double up;

        for(i=1; i<=XGRID; i++ )
                for( j=0; j<=YGRIDP1; j++ )
                        if( !boundary[i][j] )
                            base[i][j] += denud;


}


/*
 *  erfUplift:  Calculates uplift in the shape of an error function
 *      curve, with the high part of the curve toward the west. 
 *      ERFDIST gives the width, in km, over which the function
 *      ranges from 0 to ~0.998% of its maximum value. ERFOFFSET
 *      makes it possible to "slide" the position of the curve
 *      relative to the grid: if ERFOFFSET is positive, the curve
 *      is shifted toward grid west. If ERFOFFSET equals ERFDIST/2,
 *      the point of maximum slope on the curve falls at x=0,
 *      with an uplift rate at x=0 equal to half of the maximum
 *      uplift rate (i.e., half of UPRATE).
 */

erfUplift()
{
        int i, j;
        double up, scalefac=1.0, x,
        CellsToKilometers=0.001*dx;

        for(i=XGRIDP1; i>=1; i-- )
        {
                x = CellsToKilometers*i+erfoffset;
                up = erf(4.0*(x/erfdist)-2.0); 
                up = upliftRate-upliftRate*((up+1.0)/2.0);
                if( ngrabens ) scalefac=getGrabenFactor(i);
                for( j=0; j<=YGRIDP1; j++ )
                        base[i][j] += up*scalefac;
        }

}


erfSubsidence( tcurrent )
{
        int i, j, x;
        double sub,subx,erfscl,t,delt;

        delt = Dt/1000000.0;
        t = (tcurrent-tsubstart)*delt;
        sub = upliftRate*sqrt(t+delt)-upliftRate*sqrt(t);

        for(i=XGRIDP1; i>=1; i-- )
        {
                x = i+erfoffset;
                erfscl = erf(4.0*(x/erfdist)-2.0); 
                subx = sub-sub*((erfscl+1.0)/2.0);
                for( j=0; j<=YGRIDP1; j++ )
                        base[i][j] -= subx;
        }

}


halfdomeUplift()
{
        int i, j, x;
        double up;

        for(i=XGRIDP1; i>=0; i-- )
        {
                x = i+erfoffset;
                up = erf(4.0*(x/erfdist)-2.0); 
                up = upliftRate-upliftRate*((up+1.0)/2.0);
                for( j=0; j<=YGRIDP1; j++ )
                        base[i][j] += up;
        }

}


/*
 *  foldUplift: Simulates uplift of sinusoidal folds. Fold geometry and
 *      timing of uplift are determined randomly. At present, parameters
 *      are hard-coded. Designed for 3 folds on a 40x80 grid. Added 8/94.
 */ 

void foldUplift( t, action, seedValue, ifp )
int t, action, seedValue;
FILE *ifp;  /* ->input file (for INIT action only ) */
{
    int i, j, k, x, y;
    static int lamx[4], lamy[4], cx[4], cy[4], tstrt[4], tstop[4];
    static float ax[XGRID][4], ay[YGRID][4], uprt[4];
    float intrvl=333333.0,
        sigma=100000.0;
    float amp;  /* Fold amplitude */

    if( action==INIT )
    {
        /* Randomly determine wavelengths, position and amplitude for each 
         * of 3 folds.
         * For each one, a sine wave is computed to describe the fold
         * geometry in x and y.
         * Onset of folding is also determined randomly; duration depends on
         * fold amplitude and uplift rate.
         */
        printf( "Random fold parameters:\n");
        srand48( seedValue );
        for( i=0; i<4; i++ )
        {
            lamx[i] = ReadInt( ifp, "LAMX", 4 );
            if( lamx[i]==0 ) lamx[i] = 7.0+(int)(8.0*drand48());
            lamy[i] = ReadInt( ifp, "LAMY", 4 );
            if( lamy[i]==0 ) lamy[i] = 20.0+(int)(60.0*drand48());
            cx[i] = ReadInt( ifp, "CX", 2 );
            if( cx[i]==0 ) cx[i] = 2+10*i+(int)(6*drand48());
            cy[i] = ReadInt( ifp, "CY", 2 );
            if( cy[i]==0 ) cy[i] = 1+(int)(79*drand48());
            amp = ReadFloat( ifp, "AMP", 3 );
            if( amp<=0.0 ) amp = 1000.0 + 4000.0*drand48();
            uprt[i] = Dt*ReadFloat( ifp, "FUP", 3 );
            if( uprt[i]<=0.0 ) uprt[i] = Dt*(0.001 + 0.01*drand48());
            tstrt[i] = ReadInt( ifp, "TSTRT", 5 );
            if( tstrt[i]==0 ) 
                tstrt[i]=(int)(((0.5*intrvl-sigma)+(3-i)*intrvl+2.0*sigma*drand48())/(float)Dt);
            tstop[i] = tstrt[i] + (int)(amp/uprt[i]);
            printf( "lx %d, ly %d, cx %d, cy %d\n",lamx[i],lamy[i],cx[i],cy[i]);
            printf( "Amp: %.0f\tStart: %d\tStop: %d\n",amp,tstrt[i],tstop[i]);
            printf( "Uplift rate: %f m/yr\n", uprt[i]/Dt );
            for( j=0; j<=lamx[i]; j++ )
                ax[j][i]=0.5*(sin((TWOPI*j)/(float)lamx[i]+1.5*PI)+1.0);
            for( j=0; j<=lamy[i]; j++ )
                ay[j][i]=0.5*(sin((TWOPI*j)/(float)lamy[i]+1.5*PI)+1.0);
        }
    }
    
    else
      {
        for( k=0; k<4; k++ ) if( t>=tstrt[k] && t<tstop[k] ) 
            for( i=0; i<=lamx[k]; i++ )
                for( j=0; j<=lamy[k]; j++ )
                  {
                    x=i+cx[k]-(int)(0.5*lamx[k]);
                    y=j+cy[k]-(int)(0.5*lamy[k]);
                    if( y>YGRID ) y=y-YGRID;
                    else if( y<1 ) y=y+YGRID;
                    if( x>0 && x<=XGRID ) base[x][y]+=uprt[k]*ax[i][k]*ay[j][k];
                  }

        /* Compute background regional uplift */
        blockUplift( t );
      }



}


/*
**  StrikeSlip
**
**  StrikeSlip mimics strike-slip fault motion by shifting points on the grid
**  equal to or less than _fault_ horizontally. Essentially, the function
**  computes the cumulative offset since the last "shift", and if the offset
**  is greater than or equal to half a pixel, the grid points are shifted to
**  the right (toward higher Y coordinates) with respect to the others.
**  Written 10/97.
**
*/
void StrikeSlip()
{
  static float horizOffset = 0.0;
  float slipAmount = sliprate*(Dt/dx); /* Slip amt in cells per t-step*/
  int i,j;

  /* Update cumulative offset */
  horizOffset += slipAmount;

  /* If >=0.5, shift cells to mimic a 1-cell offset */
  if( horizOffset >= 0.5 )
    {
      horizOffset = horizOffset - 1.0;
      for( j=YGRIDP1; j>=2; j-- )
          for( i=1; i<=fault; i++ )
	    {
              elev[i][j] = elev[i][j-1];
              rock[i][j] = rock[i][j-1];
              chansed[i][j] = chansed[i][j-1];
              rockType[i][j] = rockType[i][j-1];
              base[i][j] = base[i][j-1];
	    }
      for( i=1; i<=fault; i++ )
	{
	  elev[i][1] = elev[i][YGRIDP1];
	  rock[i][1] = rock[i][YGRIDP1];
          chansed[i][1] = chansed[i][YGRIDP1];
          rockType[i][1] = rockType[i][YGRIDP1];
	  base[i][1] = base[i][YGRIDP1];
	  elev[i][0] = elev[i][1];
	  rock[i][0] = rock[i][1];
          chansed[i][0] = chansed[i][1];
          rockType[i][0] = rockType[i][1];
	  base[i][0] = base[i][1];
	  elev[i][YGRIDP1] = elev[i][YGRID];
	  rock[i][YGRIDP1] = rock[i][YGRID];
          chansed[i][YGRIDP1] = chansed[i][YGRID];
          rockType[i][YGRIDP1] = rockType[i][YGRID];
	  base[i][YGRIDP1] = base[i][YGRID];
	}
      streamCaptureOccurred = TRUE;
    }

    plateauUplift();

}



void Flex()
{
        int i, j, k;
        double thickl[XGRID+1][YGRIDP1], tot_thickl, tot_thicks, 
                east_edge, w0, wn, x;

        printf("Flex()...");
        
        /*
         *  Calculate loads.
         *  This assumes uniform density in rock layers.
         */
        w[0] = 0.0;
        for( i=1; i<=XGRID; i++ )
        {
                w[i] = 0.0;
                tot_thickl = 0.0;       
                tot_thicks = 0.0;       
                for( j=1; j<=YGRID; j++ )
                {
                        thickl[i][j] = rock[i][j];
                        for( k=n_layers; k>rockType[i][j]; k-- )
                                thickl[i][j] += layer[k].T;
                        tot_thickl += thickl[i][j];     
                        tot_thicks += chansed[i][j];
                }
                
                /*
                 *  Add density and gravity terms, and average in N-S
                 */
                V[i] = 
((tot_thickl*rhoc+tot_thicks*rhos)*GRAV*dx)/(double)YGRID
;
        }

        V[XGRIDP1]=V[XGRID];

        /*
         *  For floating east boundary, set the thickness of the interior
         *  plateau equal to the average thickness of rock and sediment at
         *  the easternmost column.
         */
        if( east_boundary==NO_FLUX )
        {
                east_edge = (tot_thickl+tot_thicks)/(double)YGRID;
                if( east_edge < plat_thick )
                        plat_thick = east_edge;
        }
        
        /*
         *  Shelf sediment load gets incremented by sediment yield from
         *  the previous time step.
         */
        shelf_load += sed_yield*shelf_param*0.00000001;

        /*
         *  Calculate deflection due to load in modeled region.
         */
        for( i=1; i<=XGRIDP1; i++ ) {
                w0 = unitload*V[i]*0.00000001;
                for( j=0; j<=XGRID; j++ ) {
                        x = abs(j-i) * dx;
                        wn = w0*exp(-x/alpha)*(cos(x/alpha)+sin(x/alpha));
                        w[j] += wn;
                }
        }

        /*
         *  Add deflection from shelf and interior plateau, 
         *  and subtract initial deflection profile.
         */
        for( j=0; j<=XGRIDP1; j++ )
                w[j] += shelf_w[j]*shelf_load+plat_w[j]*plat_thick-
(w_init[j]+plat_init[j]);

        w[XGRIDP1] = w[XGRID];

        printf("done\n");
        
}



double *dvector( nl, nh )
int nl, nh;
{
    double *v;

    v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
    if(!v) problem( "Unable to allocate memory in dvector()" );
    return v-nl;
}



void problem( prob )
char *prob;
{
        FILE *pfp;

        pfp = fopen( "LemErrors", "w" );
        
        fprintf( stderr, "\nProblem: %s\n\n", prob );
        fprintf( pfp, "\nProblem: %s\n\n", prob );
        fclose( pfp );
        exit( 0 );

}


outElevData( t )
int t;
{
        int i, j;
        
        fprintf( elevfp, " %d\n", t );
        for( j=0; j<=YGRIDP1; j++ )
                for( i=0; i<=XGRIDP1; i++ )
                        fprintf( elevfp, "%.2f\n", elev[i][j] );

}

outMeanElevData( t )
int t;
{
        int i, j;
        
        fprintf( mnelevfp, " %d\n", t );
        for( j=0; j<=YGRIDP1; j++ )
                for( i=0; i<=XGRIDP1; i++ )
                        fprintf( mnelevfp, "%.2f\n", meanElev[i][j] );

}



outNetData( t )
int t;
{
        int i, j, ch, ny;
        
        fprintf( netfp, " %d\n", t );
        for( j=0; j<=YGRIDP1; j++ )
                for( i=0; i<=XGRIDP1; i++ ) {
                        if( j==1 && nbry[i][j]==YGRID )
                                ny = 0;
                        else if( j==YGRID && nbry[i][j]==1 )
                                ny = YGRIDP1;
                        else ny = nbry[i][j];
                        ch=(flood[i][j]) ? 4 : chan[i][j];
                        fprintf( netfp, "%d %d %d\n", ch,nbrx[i][j],ny);
                }       

} /* outNetData */

outChanData( t )
int t;
{
        int i, j;
        
        fprintf( chanfp, " %d\n", t );
        for( j=0; j<=YGRIDP1; j++ )
                for( i=0; i<=XGRIDP1; i++ )
                        fprintf( chanfp, "%d\n", chan[i][j]);


} /* outChanData */


outSlpData( t )
int t;
{
        int i, j;
        
        fprintf( slpfp, " %d\n", t );
        for( j=1; j<=YGRID; j++ )
                for( i=1; i<=XGRID; i++ )
                        if( slp_fail[i][j] )
                                fprintf( slpfp, "%d %d %d\n", i,j,slp_fail[i][j]
);
        fprintf( slpfp, " 0\n" );


} /* outSlpData */


outSoilData( t )
int t;
{
        int i, j;
        
        fprintf( solfp, " %d\n", t );
        for( j=0; j<=YGRIDP1; j++ )
                for( i=0; i<=XGRIDP1; i++ )
                        fprintf( solfp, "%.2f\n", (chansed[i][j]+sla[i][j]+slb[i][j])/3.0 );

} /* outSoilData */


outQData( t )
int t;
{
    int i, j;
        
    fprintf( afp, " %d\n", t );
    for( j=0; j<=YGRIDP1; j++ )
        for( i=0; i<=XGRIDP1; i++ )
            fprintf( afp, "%d\n", area[i][j] );

    if( optSat )
      {
        fprintf( qfp, " %d\n", t );
        for( j=0; j<=YGRIDP1; j++ )
                for( i=0; i<=XGRIDP1; i++ )
                        fprintf( qfp, "%f\n", disch[i][j] );
      }
} 


outLithologyData( t )
int t;
{
        int i, j;
        
        fprintf( lithfp, " %d\n", t );
        for( j=1; j<=YGRID; j++ )
                for( i=1; i<=XGRID; i++ )
                        fprintf( lithfp, "%d\n", rockType[i][j] );

} /* outLithologyData */


outLayer1Data( t )
int t;
{
        int i, j;
        
        fprintf( lyrfp, " %d\n", t );
        for( j=0; j<=YGRIDP1; j++ )
                for( i=0; i<=XGRIDP1; i++ )
                        fprintf( lyrfp, "%.1f\n", rock[i][j] );

} /* outLayer1Data */


outSedData( t )
int t;
{
    int i, j;

    fprintf( sedfp, "%d\n", t );
    for( j=1; j<=YGRID; j++ )
    {
        fprintf( sedfp, "%.2f\n", (load[0][j]*cellarea)/dt_min );
        fprintf( sedfp, "%.2f\n", (load[XGRIDP1][j]*cellarea)/dt_min );
    }
}


outDelzData( t )
int t;
{
        int i, j;

        fprintf( denudfp, " %d\n", t );
        for( j=0; j<=YGRIDP1; j++ )
                for( i=0; i<=XGRIDP1; i++ )
                        fprintf( denudfp, "%f\n", (delz[i][j][0]+delz[i][j][1]
+delz[i][j][2])/(3.0*Dt) );

} /* outDenudData */


outTauData( t )
int t;
{
        int i, j;

        fprintf( taufp, " %d\n", t );
        for( j=0; j<=YGRIDP1; j++ )
                for( i=0; i<=XGRIDP1; i++ )
                        fprintf( taufp, "%.4f\n", tau[i][j] );

} /* outTauData */


outFlexData( t )
int t;
{
        int i;

        fprintf( flexfp, " %d\n", t );
        for( i=0; i<=XGRIDP1; i++ )
                fprintf( flexfp, "%.2f\n", -w[i] );
        if( t > 0 )
                for( i=0; i<=XGRIDP1; i++ )
                        fprintf( flexfp, "%.2f\n", -shelf_load*shelf_w[i] );
        else
                for( i=0; i<=XGRIDP1; i++ )
                        fprintf( flexfp, "%.2f\n", w_init[i] );
        for( i=0; i<=XGRIDP1; i++ )
                fprintf( flexfp, "%.2f\n", plat_init[i]-plat_w[i]*plat_thick );
}


char *date()
{
        char *str, *ctime();
        long time(), nseconds;

        nseconds = time((long*)0);
        str = ctime(&nseconds);
        return( str );
}


bugdump( x, y )
int x;
int y;
{
        int nx, ny, i;

        nx = nbrx[x][y];
        ny = nbry[x][y];
        printf( "Endless loop in Stream Trace: look for %d,%d\n", x, y );
        outNetData( -1 );
    for( i=1; i<=10; i++ )
    {
        printf("-> %d %d    h=%f   s=%f\n",x,y,elev[x][y],slope[x][y]);
        nx = nbrx[x][y];
        ny = nbry[x][y];
        x = nx;
        y = ny;
    }

        exit( 0 );

} /* bugdump */


void WriteRestartFile( t )
int t;
{
        FILE *fp;
        int i,j;

        if( (fp=fopen(restart_name, "w"))==NULL )
                problem( "Unable to write restart file.\n");

        fprintf( fp, "%d\n", t );
        for( j=0; j<=YGRIDP1; j++ )
                for( i=0; i<=XGRIDP1; i++ )
                        fprintf( fp, "%d\t%.2f\t%.2f\t%.2f\n",
                                rockType[i][j],rock[i][j],base[i][j],chansed[i][j])
;
        if( flexure )
        {
                fprintf( fp, "%.2f\n", plat_thick );
                fprintf( fp, "%.2f\n", shelf_load );
                for( i=0; i<=XGRIDP1; i++ )
                        fprintf( fp, "%.2f\n", w_init[i] );
        }

/* For unknown reasons this fclose statement was producing segmentation
faults, during testing of variable boundary version. Began happening
after a variable tmpCharArray[][] was added to ReadLithology(). 
        fclose( fp );*/
                

}


/* dmatrix: Allocates a double matrix with range [nrl..nrh][ncl..nch] */

double **dmatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
    int i;
    double **m;

    m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
    if(!m) problem("allocation failure 1 in dmatrix()" );
    m -= nrl;

    /* Allocate rows and set pointers to them */
    for(i=nrl;i<=nrh;i++) {
        m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
        if(!m[i]) problem("allocation failure 2 in dmatrix()");
        m[i] -= ncl;
    }

    return m;

} /* dmatrix */


int **imatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
    int i, **m;

    m=(int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
    if(!m) problem("allocation failure 1 in imatrix()");
    m -= nrl;
    for(i=nrl;i<=nrh;i++) {
        m[i]=(int *)malloc((unsigned) (nch-ncl+1)*sizeof(int));
        if(!m[i]) problem("allocation failure 2 in imatrix()");
        m[i] -= ncl;
    }
    return m;
}


char **cmatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
    int i; 
    char **m;

    m=(char **)malloc((unsigned) (nrh-nrl+1)*sizeof(char*));
    if(!m) problem("allocation failure 1 in cmatrix()");
    m -= nrl;
    for(i=nrl;i<=nrh;i++) {
        m[i]=(char *)malloc((unsigned) (nch-ncl+1)*sizeof(char));
        if(!m[i]) problem("allocation failure 2 in imatrix()");
        m[i] -= ncl;
    }
    return m;
}


/* free_dmatrix: Frees a matrix allocated with dmatrix() */

void free_dmatrix(m,nrl,nrh,ncl,nch)
double **m;
int nrl,nrh,ncl,nch;
{
    int i;

    for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
    free((char*) (m+nrl));
}


void free_imatrix(m,nrl,nrh,ncl,nch)
int **m;
int nrl,nrh,ncl,nch;
{
    int i;

    for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
    free((char*) (m+nrl));
}

/*** End of GOLEM.C ***/
