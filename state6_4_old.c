/*
 * State elements located at subsampled locations along
 * boundaries, always integer subdivision.  State 
 * definitions can represent local eigenvectors to
 * improve conditioning.
 *
 * Paul Fieguth
 * July, 2003
 *
 * Bugs:
 * *
 * Modifications:
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifndef Matlab5
#include "cmex.h"
#else
#include "mex.h"
#endif
#include "smth.h"
#include "clock.h"

#include "statecorr.h"

#ifndef Simple_State
#error Variable Simple_State must be defined.
#endif

#ifndef Matlab_Callback
#error Matlab Callbacks must be defined for State6_4.c
#endif

#ifndef CorrFn_Region
#error Correlation function regions must be specified (CorrFn_Region)
#endif

#define MaxStateSize 2800
#define MaxBlockSize 16
#define MaxEigEls    6
#define MaxEigs      16

typedef struct {
  int len;
  int *dx;
  int *dy;
  double *statevals;
} state;

typedef struct {
  int f, sc;
  int x, y;
  int wx, wy;
  int dirn;
  int *posx;
  int *posy;
  state **sp;
  int *reg;
} sinfo;

static int num_fields;
static int xstate, ystate, fstate, xsiz, ysiz;

static int st_samp[Max_Fields][Max_Scales][4];
static int st_sizx[Max_Fields][Max_Scales][4];
static int st_sizy[Max_Fields][Max_Scales][4];
static int st_els[Max_Fields][Max_Scales][4];
static double st_sub[Max_Fields][Max_Scales][4];
static int st_eig[Max_Fields][Max_Scales][4];
static double st_elim[Max_Fields][Max_Scales][4];

static state *statedef[MaxBlockSize+1][MaxBlockSize+1][MaxEigEls+1][MaxEigs+1];

/* local, custom state storage */
#define NumTempState (5*Max_Fields*MaxEigs)
static state temp_states[NumTempState];
static int temp_ind;

/* parameters available from ms_cc6_corr */
extern int scales;
extern int ovlp[];
extern int idx[], idy[];
extern int debug_mode, debug_opt;
extern int clk_ovlp, clk_dim, clk_cov;
extern int realloc_memory;

/*
 * Sampled Boundary routine
 */
int mexStateFn(nrhs, prhs)
     int    nrhs;
     Matrix *prhs[];
{
  int i, j, k, l, f, rowofs;

  /* Check for proper number of arguments */
  num_fields = 0;
  if (nrhs > 0) num_fields = ((int) floor(mxGetScalar( prhs[0] )+0.5));
  if (nrhs < 1+num_fields) {
    mexPrintf( "\nCorrect function usage: (file: state6_4.c)\n" );
    mexPrintf( "   State Parameters:  fields, prm(1), ..., prm(fields)]\n" );
    mexPrintf( "    fields - number of separate random fields\n" );
    mexPrintf( "             fields(2) - Horz boundary:  0 asserted, 1 periodic, 2 empty\n" );
    mexPrintf( "             fields(3) - Vert boundary:  0 asserted, 1 periodic, 2 empty\n" );
    mexPrintf( "             fields(4) - If present, asserts independent fields\n" );
    mexPrintf( "    prmi   - state boundary parameters for field i\n" );
    mexPrintf( "             scalar - boundary density (integer) \n" );
    mexPrintf( "             row    - boundary density by scale (coarse to fine)\n" );
    mexPrintf( "             5rows  - density (# samples per side length, negative for global\n" );
    mexPrintf( "                      size x, size y (integer)\n" );
    mexPrintf( "                      a.b:  a - num samp elements (sparsification of eigenvector)\n" );
    mexPrintf( "                            b - fraction of pix to elim (sparsif. of domain)\n" );
    mexPrintf( "                      a.b:  a - num eig (integer)\n" );
    mexPrintf( "                            b - relative eval to keep (zero - keep all\n" );
    mexPrintf( "            10rows  - first 5 rows - horizontal direction\n" );
    mexPrintf( "                    - next  5 rows - vertical direction\n" );
    mexPrintf( "            20rows  - 5 horz, 5 vert, 5 corners, 5 joints\n" );
    
    if (prhs != NULL) mexErrMsgTxt( "Wrong Number of Parameters Given." );
    return(0);
  }

  if (num_fields > Max_Fields)
    mexErrMsgTxt( "Too many fields requested." );

  xstate = ystate = fstate = 0;
  if (mxGetM(prhs[0])+mxGetN(prhs[0])>2)
    xstate = ((int) floor((mxGetPr(prhs[0]))[1]+0.5));
  if (mxGetM(prhs[0])+mxGetN(prhs[0])>3)
    ystate = ((int) floor((mxGetPr(prhs[0]))[2]+0.5));
  if (mxGetM(prhs[0])+mxGetN(prhs[0])>4)
    fstate = 1;

  for (f=0; f<num_fields; f++) {
    for (i=0; i<Max_Scales; i++) 
      for (j=0; j<4; j++) {
	st_samp[f][i][j] = floor((mxGetPr(prhs[1+f]))[min(i,mxGetN(prhs[1+f])-1)*mxGetM(prhs[1+f])]+0.5);
	st_sizx[f][i][j] = 1;
	st_sizy[f][i][j] = 1;
	st_els[f][i][j] = 0;
	st_sub[f][i][j] = 0;
	st_eig[f][i][j] = 1;
	st_elim[f][i][j] = 0;
      }

    if (mxGetM(prhs[1+f])>4) {
      for (j=0; j<4; j++) {
	if (mxGetM(prhs[1+f])>j*5+4) rowofs = j*5;

	for (i=0; i<mxGetN(prhs[1+f]); i++) {
	  st_samp[f][i][j] = floor((mxGetPr(prhs[1+f]))[rowofs+i*mxGetM(prhs[1+f])]+0.5);
	  st_sizx[f][i][j] = floor((mxGetPr(prhs[1+f]))[rowofs+1+i*mxGetM(prhs[1+f])]+0.5);
	  st_sizy[f][i][j] = floor((mxGetPr(prhs[1+f]))[rowofs+2+i*mxGetM(prhs[1+f])]+0.5);
	  st_els[f][i][j] = floor((mxGetPr(prhs[1+f]))[rowofs+3+i*mxGetM(prhs[1+f])]+0.0001);
	  st_sub[f][i][j] = (mxGetPr(prhs[1+f]))[rowofs+3+i*mxGetM(prhs[1+f])]-st_els[f][i][j];
	  st_eig[f][i][j] = floor((mxGetPr(prhs[1+f]))[rowofs+4+i*mxGetM(prhs[1+f])]+0.0001);
	  st_elim[f][i][j] = (mxGetPr(prhs[1+f]))[rowofs+4+i*mxGetM(prhs[1+f])]-st_eig[f][i][j];
	}
      }
    }

    /* assert limits */
    for (i=0; i<Max_Scales; i++) {
      if ((st_samp[f][i][0] + 0.5) * (st_samp[f][i][1] + 0.5) < 0)
	mexErrMsgTxt( "Cannot mix local and global modes in vert/horz sampling." );
      if (st_samp[f][i][0] > -0.5) {
	for (j=0; j<4; j++) {
	  /* size must be odd */
	  if ((st_sizx[f][i][j] & 1) == 0) st_sizx[f][i][j]--;
	  st_sizx[f][i][j] = max(1,st_sizx[f][i][j]);

	  if ((st_sizy[f][i][j] & 1) == 0) st_sizy[f][i][j]--;
	  st_sizy[f][i][j] = max(1,st_sizy[f][i][j]);
	}
      }

      for (j=0; j<4; j++) {
	st_els[f][i][j] = max(0,min(st_els[f][i][j],MaxEigEls-1));
	st_eig[f][i][j] = max(1,min(st_eig[f][i][j],MaxEigs-1));
      }
    }

    /* assert finest scales */
    st_samp[f][scales-2][0] = -1;
    st_samp[f][scales-1][0] = -1;
  }

  /* initialize state definitions */
  if (realloc_memory) {
    for (i=0;i<MaxBlockSize+1;i++) for (j=0;j<(MaxBlockSize+1);j++)
      for (k=0;k<(MaxEigEls+1);k++) for (l=0;l<(MaxEigs+1);l++)
	statedef[i][j][k][l] = NULL;
  }
  state_gen( 1, 1, 0, 0, 0, 0 );

  /* clear temp state memory */
  for (temp_ind=0; temp_ind<NumTempState; temp_ind++) {
    temp_states[temp_ind].len = 0;
    temp_states[temp_ind].dx = NULL;
    temp_states[temp_ind].dy = NULL;
    temp_states[temp_ind].statevals = NULL;
  }
  temp_ind = 0;

  /* get overall domain size */
  overlap_data( 0, 0, scales, ovlp, idx, idy, &i, &i, &xsiz, &ysiz );

  return( 1+num_fields );
}


/* check that all of defined state element lies within single region */
int state_check( int f, int x, int y, state *s )
{
  int i;

  for (i=1;i<s->len;i++)
    if (!(corr_pair(f,x+s->dx[0],y+s->dy[0],x+s->dx[i],y+s->dy[i])))
      return(0);
  return(1);
}

void state_gen( int wx, int wy, int numels, int numeig, double evallim, double sparsity )
{
  char matlabstr[500];
  double *px, *py, *ps, *pl, d_zero, d_one;
  int ind, i, j, l;

  d_zero = 0; d_one = 1;
  if ((wx == 1) && (wy == 1)) {
    px = &d_zero;
    py = &d_zero;
    pl = &d_one;
    ps = &d_one;
    numeig = 1;
  }
  else {
    mexEvalString( "st64 = exist('state_gen','file');" );
    if (mxGetScalar( mexGetArrayPtr( "st64", "caller" ) ) != 2)
      mexErrMsgTxt( "Unable to find Matlab support function state_gen.m" );
    
    sprintf( matlabstr, "[st64x,st64y,st64l,st64s] = state_gen(%d,%d,%d,%d,%lf,%lf);", wx, wy, numels, numeig, evallim, sparsity );
    mexEvalString( matlabstr );
  
    px = mxGetPr( mexGetArrayPtr( "st64x", "caller" ) );
    py = mxGetPr( mexGetArrayPtr( "st64y", "caller" ) );
    pl = mxGetPr( mexGetArrayPtr( "st64l", "caller" ) );
    ps = mxGetPr( mexGetArrayPtr( "st64s", "caller" ) );
  }

  for (ind=0,i=0;i<numeig;i++) {
    if ((l = (int)pl[i]) > 0) {
      /* have a state - allocate and copy */
      statedef[wx][wy][numels][i] = pf_mxCalloc( 1, sizeof(state) );
      statedef[wx][wy][numels][i]->len = l;
      statedef[wx][wy][numels][i]->dx = pf_mxCalloc( l, sizeof(int) );
      statedef[wx][wy][numels][i]->dy = pf_mxCalloc( l, sizeof(int) );
      statedef[wx][wy][numels][i]->statevals = pf_mxCalloc( l, sizeof(double) );

      for (j=0;j<l;j++) {
	statedef[wx][wy][numels][i]->dx[j] = (int) px[ind];
	statedef[wx][wy][numels][i]->dy[j] = (int) py[ind];
	statedef[wx][wy][numels][i]->statevals[j] = ps[ind++];
      }
    }
  }
}

/* interpolate larger state from smaller one */
state *state_interp( int f, int sc, int wx, int wy, int swx, int swy, int eig, int type )
{
  state *sint, *stmp;
  int numels, l;

  /* generate state if needed */
  numels = st_els[f][sc][type];
  if (statedef[swx][swy][numels][0] == NULL) {
    state_gen(swx,swy,numels,st_eig[f][sc][type],st_elim[f][sc][type],st_sub[f][sc][type]);
  }

  sint = statedef[swx][swy][numels][eig];
  if (sint == NULL) return(NULL);
  stmp = temp_states+temp_ind;

  /* free/allocate as needed */
  if (stmp->len < sint->len) {
    if (stmp->dx != NULL) mxFree( stmp->dx );
    if (stmp->dy != NULL) mxFree( stmp->dy );
    if (stmp->statevals != NULL) mxFree( stmp->statevals );

    stmp->dx = (int *) mxCalloc( sint->len, sizeof(int) );
    stmp->dy = (int *) mxCalloc( sint->len, sizeof(int) );
    stmp->statevals = (double *) mxCalloc( sint->len, sizeof(double) );
  }
  stmp->len = sint->len;

  /* interpolate over non-zero terms */
  for (l=0; l<sint->len; l++) {
    stmp->dx[l] = (sint->dx[l] * wx)/swx;
    stmp->dy[l] = (sint->dy[l] * wy)/swy;
    stmp->statevals[l] = sint->statevals[l];
  }

  temp_ind = (temp_ind + 1) % NumTempState;
  return(stmp);
}

int state_glob( int f, int sc, int x, int y, int wx, int wy, int type,
		int regset, int *posx, int *posy, state **sp, int *reg )
{
  int i, swx, swy, acteig, numeig, numels;

  numeig = st_eig[f][sc][type];
  numels = st_els[f][sc][type];

  if ((wx <= MaxBlockSize) && (wy <= MaxBlockSize)) {
    if (statedef[wx][wy][numels][0] == NULL)
      state_gen(wx,wy,numels,st_eig[f][sc][type],st_elim[f][sc][type],st_sub[f][sc][type]);

    for (acteig=0,i=0;i<numeig;i++) 
      if (statedef[wx][wy][numels][i] != NULL) {
	posx[acteig] = x; posy[acteig] = y; reg[acteig] = regset;
	sp[acteig] = statedef[wx][wy][numels][i];
	acteig++;
      }
    return(acteig);
  }
      
  /* will need to interpolate out from smaller state */
  swx = min(MaxBlockSize,wx);
  swy = min(MaxBlockSize,wy);
  for (acteig=0,i=0;i<numeig;i++) {
    if ((sp[acteig] = state_interp( f, sc, wx, wy, swx, swy, i, type )) != NULL) {
      posx[acteig] = x; posy[acteig] = y; reg[acteig] = regset;
      acteig++;
    }
  }
  return(acteig);
}

int state_group( int f, int sc, int x, int y, int wx, int wy, int type, 
		 int *posx, int *posy, state **sp, int *reg, int nonull )
{
  int g, i, j, ind, lind, cnt;
  int corr_region( int f, int x, int y, int wx, int wy ); /* copied from 570 in corr_land.c*/
  if ((g = corr_region( f,x,y,wx,wy )) > 0)
    return( state_glob( f,sc,x,y,wx,wy,type, g, posx, posy, sp, reg ) );

  if (g==0) {
    if (nonull) return(0);
    posx[0] = x; posy[0] = y; reg[0] = 0; sp[0] = statedef[1][1][0][0]; 
    return(1);
  }

  ind = 0;
  for (i=0;i<wx;i++) for (j=0;j<wy;j++)
    if (!(corr_null(f,x+i,y+j))) {
      cnt = 1;
      for (lind=0;lind<ind;lind+=1)
	if (corr_types(f,reg[lind],corr_type(f,x+i,y+j)))
	  { cnt = 0; break; }
      if (cnt == 1) {
	ind += state_glob( f,sc,x,y,wx,wy,type, corr_type( f,x+i,y+j ), posx+ind, posy+ind, sp+ind, reg+ind );
      }
    }
  return( ind );
}

int state_side( int f, int sc, int x, int y, int type, int len, int elim, int *posx, int *posy, state **sp, int *reg, int nonull )
{
  int eigjoint, i, ind, s, dirn, wx, wy;

  ind = 0; 
  dirn = type;
  s = st_samp[f][sc][dirn]+1;

  for (i=1; i<s; i++) {
    type = dirn;
    if (((abs(elim)*i) % s)==0) type = 3;
    wx = st_sizx[f][sc][type];
    wy = st_sizy[f][sc][type];
    if ((elim < 0) || (type < 2))
      ind += state_group( f, sc, x+(1-dirn)*(len*i)/s-(wx-1)/2, y+dirn*(len*i)/s-(wy-1)/2, wx, wy, type, posx+ind, posy+ind, sp+ind, reg+ind, nonull+ind );
  }
  
  return( ind );
}

/* horz segments 0, 1, 4;  vert segments 2, 3, 5 */
static int dirn[6] = { 0, 0, 1, 1, 0, 1 };

int state_dim( int f, int sc, int x, int y, int wx, int wy, int dx, int dy, 
	       int *posx, int *posy, state **sp, int *reg )
{
  int i, j, ind, cnt;
  int cwx, cwxo, cwy, cwyo;

  /* handle finest scale explicitly */
  if ((wx==1) && (wy==1)) { 
    posx[0] = x; posy[0] = y; sp[0] = statedef[1][1][0][0]; reg[0] = corr_type(f,x,y);
    return(1);
  }

  if (st_samp[f][sc][0] < 0) 
    /* global, full representation of state, rather than boundary */
    return( state_group( f, sc, x, y, wx, wy, 0, posx, posy, sp, reg, 0 ) );

  /* have regular boundary, spaced around outside and inside of domain */
  ind = 0;
  cwx = st_sizx[f][sc][2]; cwxo = -(cwx-1)/2;
  cwy = st_sizy[f][sc][2]; cwyo = -(cwy-1)/2;

  /* always need top corner unless boundary empty */
  if (!((x == 0) && (y == 0) && (xstate == ST_EMPTY) && (ystate == ST_EMPTY)))
    ind += state_group( f, sc, x+cwxo, y+cwyo, cwx, cwy, 2, posx+ind, posy+ind, sp+ind, reg+ind, ind );

  /* want left / top sides unless boundary empty */
  if (!((y == 0) && (ystate == ST_EMPTY)))
    ind += state_side( f, sc, x, y, 0, wx, -dy, posx+ind, posy+ind, sp+ind, reg+ind, ind );
  if (!((x == 0) && (xstate == ST_EMPTY)))
    ind += state_side( f, sc, x, y, 1, wy, -dx, posx+ind, posy+ind, sp+ind, reg+ind, ind );

  /* other states needed if internal or regular boundary */
  if (((x+wx<xsiz) || (xstate == ST_NORM)) && ((y>0) || (ystate != ST_EMPTY)))
    ind += state_group( f, sc, x+wx-1+cwxo, y+cwyo, cwx, cwy, 2, posx+ind, posy+ind, sp+ind, reg+ind, ind );
  if ((x+wx<xsiz) || (xstate == ST_NORM))
    ind += state_side( f, sc, x+wx-1, y, 1, wy, -dx, posx+ind, posy+ind, sp+ind, reg+ind, ind  );

  if (((y+wy<ysiz) || (ystate == ST_NORM)) && ((x>0) || (xstate != ST_EMPTY)))
    ind += state_group( f, sc, x+cwxo, y+wy-1+cwyo, cwx, cwy, 2, posx+ind, posy+ind, sp+ind, reg+ind, ind );
  if ((y+wy<ysiz) || (ystate == ST_NORM))
    ind += state_side( f, sc, x, y+wy-1, 0, wx, -dy, posx+ind, posy+ind, sp+ind, reg+ind, ind  );

  if (((x+wx<xsiz) || (xstate == ST_NORM)) && ((y+wy<ysiz) || (ystate == ST_EMPTY)))
    ind += state_group( f, sc, x+wx-1+cwxo, y+wy-1+cwyo, cwx, cwy, 2, posx+ind, posy+ind, sp+ind, reg+ind, ind );

  /* need internal bounaries - first horizontal */
  for (i=1; i<dy; i++)
    ind += state_side( f, sc, x, y+i*wy/dy, 0, wx, -dy, posx+ind, posy+ind, sp+ind, reg+ind, ind  );

  /* next vertical, eliminating crossing pixels if present */
  for (i=1; i<dx; i++)
    ind += state_side( f, sc, x+i*wx/dx, y, 1, wy, dx, posx+ind, posy+ind, sp+ind, reg+ind, ind  );

  return(ind);
}

/* assume state dimension never exceeds MaxStateSize elements */
static int posx[2*MaxStateSize], posy[2*MaxStateSize], fld[2*MaxStateSize];
static int region[2*MaxStateSize];
static state *states[2*MaxStateSize];

void statefn( double (*corrfn)(), 
	      int sc1, int s1, int x1, int y1, int wx1, int wy1,
	      int sc2, int s2, int x2, int y2, int wx2, int wy2,
	      int *state_dim1, int *state_dim2, double *cov ) 
{
  state *st1, *st2;
  int f, dim, d1, d2, sd1, sd2, si1, si2, covptr;
  int dx1, dx2, dy1, dy2, l1, l2, ps1, ps2;

  time_track( CLK_START, clk_dim );

  if (cov == NULL) {
    if (state_dim1 != NULL) {
      if (sc1 >= scales-1) {
	*state_dim1 = num_fields;
      }
      else {
	for (dim=0,f=0; f<num_fields; f++) 
	  dim += state_dim( f, sc1, x1, y1, wx1, wy1, idx[sc1], idy[sc1], posx, posy, states, region );
	*state_dim1 = dim;
if (dim==0) printf("Zero state %d %d %d %d %d %d\n",sc1,s1,x1,y1,wx1,wy1);
      }
    }
    if (state_dim2 != NULL) {
      if (sc2 >= scales-1) {
	*state_dim2 = num_fields;
      }
      else {
	for (dim=0,f=0; f<num_fields; f++) 
	  dim += state_dim( f, sc2, x2, y2, wx2, wy2, idx[sc2], idy[sc2], posx, posy, states, region );
	*state_dim2 = dim;
if (dim==0) printf("Zero state %d %d %d %d %d %d\n",sc2,s2,x2,y2,wx2,wy2);
      }
    }
    time_track( CLK_STOP, clk_dim );
    return;
  }

  /* need a covariance - at finest scale? */
  if ((sc1 == scales-1) && (sc2 == scales-1)) {
    time_track( CLK_STOP, clk_dim );
    time_track( CLK_START, clk_cov );
    
    for (covptr=0,d2=0; d2<num_fields; d2++)
      for (d1=0; d1<num_fields; d1++,covptr++) {
	if ((d1 != d2) && (fstate != 0)) { 
	  cov[covptr] = 0.0; 
	}
	else
	  cov[covptr] = corrfn( d1, x1, y1, 0, d2, x2, y2, 0 );
      }
    
    time_track( CLK_STOP, clk_cov );
    return;
  }

  /* compute state 1 arrangement */
  sd1 = dim = 0;
  for (f=0; f<num_fields; f++) {
    dim += state_dim( f, sc1, x1, y1, wx1, wy1, idx[sc1], idy[sc1], posx+sd1, posy+sd1, states+sd1, region+sd1 );
    while (sd1<dim) fld[sd1++] = f;
  }
  if (state_dim1 != NULL) *state_dim1 = dim;
  time_track( CLK_STOP, clk_dim );

  /* print state if desired */
  if (debug_mode && debug_opt) {
    printf( "State Region (%d,%d) to (%d,%d)\n",x1,y1,x1+wx1,y1+wy1);
    for (d1=0; d1<sd1; d1++) {
      printf( "F %d (%d,%d) region %d:  ", fld[d1],posx[d1],posy[d1],region[d1] );
      for (si1=0; si1<states[d1]->len; si1++)
	printf( "(%d,%d)%4.2lf  ",states[d1]->dx[si1],states[d1]->dy[si1],states[d1]->statevals[si1] );
      printf( "\n" );
    }
  }

  /* plot state if desired */
#ifdef Matlab_Callback
  if (debug_mode && debug_opt) {
    char temps[500];
    mexEvalString( "hold on;" );
    
    for (f=0;f<num_fields;f++) {
      sprintf( temps, "plot([%d %d %d %d %d],[%d %d %d %d %d]);", x1+(wx1+2)*f,x1+wx1+(wx1+2)*f,x1+wx1+(wx1+2)*f,x1+(wx1+2)*f,x1+(wx1+2)*f,y1,y1,y1+wy1,y1+wy1,y1);
      mexEvalString( temps );
    }

    for (d1=0; d1<sd1; d1++) {
      for (si1=0; si1<states[d1]->len; si1++) {
	sprintf( temps, "text(%d,%d,'%c');", posx[d1]+(wx1+2)*fld[d1]+states[d1]->dx[si1], posy[d1]+states[d1]->dy[si1], (char) ('0'+(d1 % 10)) );
	mexEvalString( temps );
      }
    }
  }
#endif

  
  /* if self covariance, take advantage of symmetry */
  if ((sc1==sc2) && (s1==s2)) {
    time_track( CLK_START, clk_cov );
    covptr = 0;
    for (d2=0; d2<sd1; d2++)  for (d1=0; d1<sd1; d1++) {
      if (d2>d1) { 
	/* symmetric element, don't recompute */
	cov[covptr++] = cov[d2+d1*sd1]; 
      }
      else {
	/* new element - will have to generate */
	cov[covptr] = 0.0;
	if ((d1==d2) || (((fld[d2] == fld[d1]) || (fstate == 0)) && (corr_types(fld[d1],region[d1],region[d2])))) {
	  st1 = states[d1]; st2 = states[d2];
	  for (l2=0; l2<st2->len; l2++)
	    for (l1=0; l1<st1->len; l1++)
	      cov[covptr] += st1->statevals[l1] * st2->statevals[l2] * corrfn_quick( fld[d1], posx[d1]+st1->dx[l1], posy[d1]+st1->dy[l1], posx[d2]+st2->dx[l2], posy[d2]+st2->dy[l2] );
	}
	covptr++;
      }
    }
  }
  else {
    /* compute state 2 arrangement */
    time_track( CLK_START, clk_dim );
    sd2 = dim;
    for (f=0; f<num_fields; f++) {
      dim += state_dim( f, sc2, x2, y2, wx2, wy2, idx[sc2], idy[sc2], posx+sd2, posy+sd2, states+sd2, region+sd2 );
      while (sd2<dim) fld[sd2++] = f;
    }
    if (state_dim2 != NULL) *state_dim2 = sd2-sd1;
    time_track( CLK_STOP, clk_dim );

    /* matrix is not symmetric; generate full */
    time_track( CLK_START, clk_cov );
    covptr = 0;
    for (d2=sd1; d2<sd2; d2++)  for (d1=0; d1<sd1; d1++) {
      cov[covptr] = 0.0;
      if (((fld[d2] == fld[d1]) || (fstate == 0)) && (corr_types(fld[d1],region[d1],region[d2]))) {
	st1 = states[d1]; st2 = states[d2];
	for (l2=0; l2<st2->len; l2++)
	  for (l1=0; l1<st1->len; l1++)
	    cov[covptr] += st1->statevals[l1] * st2->statevals[l2] * corrfn_quick( fld[d1], posx[d1]+st1->dx[l1], posy[d1]+st1->dy[l1], posx[d2]+st2->dx[l2], posy[d2]+st2->dy[l2] );
      }
      covptr++;
    }
  }

  time_track( CLK_STOP, clk_cov );
  return;
}
