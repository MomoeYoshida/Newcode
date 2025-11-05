/*
 * Oxford / OI correlation function.
 * Allows rescaled variances by pixel to function as a space-varying
 * prior.
 * Also allows separated bodies of water to be estimated independently
 * by decoupling, based on a matrix of specified regions.
 *
 * Compiler directive:
 *    Detailed_Distance - include matrix of local details
 *    Corrfn_Region     - allow specification of region to corrfn
 *
 * Paul Fieguth
 * October, 1998
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

#include "statecorr.h"

/* locally-declared correlation function */
#ifdef Detailed_Distance
double correl_stat(int,int,int,double);
double correl_nonstat( int f, int dx, int dy, int halfind, double d );
#else
double correl_stat(int,int,int);
double correl_nonstat( int f, int dx, int dy, int halfind );
#endif
double correl_fn( int f, double r );
double correl_length(int,double);

static int corrtype[Max_Fields];
static int wraptype[Max_Fields];
static int sigma_type[Max_Fields];
static double *sgm[Max_Fields];
static double len_x[Max_Fields], len_y[Max_Fields];
static double tp = 2*3.14159265;
static int num_fields;

static double *landptr, *distptr, *adjptr[Max_Fields];
static int xsiz, ysiz, adjsiz;

extern int debug_mode;
extern int debug_opt;
extern int scales;
extern int idx[], idy[], ovlp[];

#define LocalPix 33
double localcorrs[Max_Fields][2*LocalPix+1][2*LocalPix+1];
double localdists[2*LocalPix+1][2*LocalPix+1];
double *usercorrs[Max_Fields];
int    userorigin[Max_Fields][2];
int    usersize[Max_Fields];

/* variables for nonstationary correlation */
double *cl11[Max_Fields], *cl12[Max_Fields], *cl22[Max_Fields];
int corr_spvar[Max_Fields], corr_anis[Max_Fields];


/*
 * Exponential Correlation routine
 */
int mexCorrFn(nrhs, prhs)
     int    nrhs;
     Matrix *prhs[];
{
  int f, i, j, pofs;

  pofs = 0;
#ifdef Detailed_Distance
  pofs = 1;
#endif

  /* Check for proper number of arguments */
  num_fields = 1;
  if (nrhs > 3+pofs) { num_fields = (int) mxGetScalar( prhs[3+pofs] ); }
  if (nrhs < 4+2*num_fields+pofs) {
    mexPrintf( "\nCorrect function usage: (file: corr_land.c)\n" );
    mexPrintf( "   Corr. Param: corrtype, region, radjacent, fields, s1, ..., sfields, \n" );
    mexPrintf( "                l1, ..., lfields\n" );
    mexPrintf( "    corrtype - correlation type (scalar, or [type_1 ... type_fields]\n" );
    mexPrintf( "               0 - Gaussian\n" );
    mexPrintf( "               1 - Separable Exponential\n" );
    mexPrintf( "               2 - NonSep Exponential\n" );
    mexPrintf( "               3 - Smooth Non-Gaussian\n" );
    mexPrintf( "               4 - Specified Stationary Correlation Matrix\n" );
    mexPrintf( "                   (size odd x odd, origin in middle of array)\n" );
    mexPrintf( "               5 - First-order modified Bessel K\n" );
    mexPrintf( "               6 - Spherical\n" );
    mexPrintf( "               7 - Logistic\n" );

    mexPrintf( "      second row, optional:\n" );
    mexPrintf( "               0 - Nonwrapping boundaries\n" );
    mexPrintf( "               1 - Wrapping boundaries (may cause singularities)\n" );
    mexPrintf( "    region - integer matrix of connected regions;\n" );
    mexPrintf( "             region 0 is not estimated - independent pixels\n" );
    mexPrintf( "             do not overlap region, si, or li arrays\n" );
#ifdef Detailed_Distance
    mexPrintf( "    distances - integer matrix of local distances\n" );
    mexPrintf( "             inter-pixel distance is larger of Euclidean and this matrix\n" );
#endif
    mexPrintf( "    radjacent - adjacency matrix (symmetric); entry (i,j)\n" );
    mexPrintf( "             sets 1-adjacent or 0-independent\n" );
    mexPrintf( "             can concat matrices horiz. for field-dependent adjacency\n" );
    mexPrintf( "    fields - number of planes to be estimated\n" );
    mexPrintf( "    si     - std. dev. of process for field i;\n" );
    mexPrintf( "             if a matrix, has spatially varying variance\n" );
    mexPrintf( "    li     - correlation length (one std. dev. spatially) of field i\n" );
    mexPrintf( "             give 2-vector [lx_i ly_i] for anisotropic corr. length\n" );
    mexPrintf( "             give matrix for space-varying correl. lengths\n" );
    mexPrintf( "             give 1x3 stacked matrix for space-varying, anisotropic\n");
    mexPrintf( "                  correl. lengths (a11,a12,a22 of 2x2 matrix inverse)\n" );
    mexPrintf( "             give matrix of correlations for option 4; origin in middle\n" );
    mexPrintf( "                  matrix must be of size oddxodd; origin should be normalized\n" );

    
    if (prhs != NULL) mexErrMsgTxt( "Wrong Number of Parameters Given." );
    return(0);
  }

  if (num_fields > Max_Fields) 
    mexErrMsgTxt( "Too many fields." );

  for (i=0; i<num_fields; i++) {
    wraptype[i] = 0;
    corrtype[i] = (int) floor((mxGetPr(prhs[0]))[min(i,mxGetN(prhs[0])-1)*mxGetM(prhs[0])]+0.5);
    if (mxGetM(prhs[0])>1)
      wraptype[i] = (int) floor((mxGetPr(prhs[0]))[1+min(i,mxGetN(prhs[0])-1)*mxGetM(prhs[0])]+0.5);
  }

  overlap_data( 0, 0, scales, ovlp, idx, idy, &i, &i, &xsiz, &ysiz );

  landptr = mxGetPr(prhs[1]);
  if ((mxGetN(prhs[1])!=xsiz) || (mxGetM(prhs[1])!=ysiz))
    mexErrMsgTxt( "Wrong sized region matrix." );
#ifdef Detailed_Distance
  distptr = mxGetPr(prhs[2]);
  if ((mxGetN(prhs[2])!=xsiz) || (mxGetM(prhs[2])!=ysiz))
    mexErrMsgTxt( "Wrong sized distance detail matrix." );
#endif

  for (i=0; i<num_fields; i++)
    adjptr[i] = mxGetPr(prhs[2+pofs]);
  adjsiz = mxGetM(prhs[2+pofs]);

  if ((mxGetM(prhs[2+pofs]) != mxGetN(prhs[2+pofs])) && (mxGetM(prhs[2+pofs])*num_fields != mxGetN(prhs[2+pofs])))
    mexErrMsgTxt( "Adjacency matrix must be square (or concatenation of square matrices)." );
  if (mxGetM(prhs[2+pofs])*num_fields == mxGetN(prhs[2+pofs]))
    for (i=0; i<num_fields; i++)
      adjptr[i] = mxGetPr(prhs[2+pofs])+i*mxGetM(prhs[2+pofs])*mxGetM(prhs[2+pofs]);

  for (i=0; i<num_fields; i++) {
    sgm[i] = mxGetPr(prhs[4+i+pofs]);
    sigma_type[i] = (mxGetM(prhs[4+i+pofs])+mxGetN(prhs[4+i+pofs])>2);

    corr_spvar[i] = corr_anis[i] = 0;
    if (corrtype[i] != 4) {
      len_x[i] = len_y[i] = mxGetScalar(prhs[4+num_fields+i+pofs]);
      if (mxGetM(prhs[4+num_fields+i+pofs])+mxGetN(prhs[4+num_fields+i+pofs])>2) {
	if (mxGetM(prhs[4+num_fields+i+pofs])*mxGetN(prhs[4+num_fields+i+pofs])==2) {
	  len_y[i] = mxGetPr(prhs[4+num_fields+i+pofs])[1];
	} else if ((mxGetM(prhs[4+num_fields+i+pofs])==ysiz) && (mxGetN(prhs[4+num_fields+i+pofs])==xsiz)) {
	  corr_spvar[i] = 1;
	  cl11[i] = mxGetPr(prhs[4+num_fields+i+pofs]);
	} else if ((mxGetM(prhs[4+num_fields+i+pofs])==ysiz) && (mxGetN(prhs[4+num_fields+i+pofs])==3*xsiz)) {
	  corr_spvar[i] = 1;
	  corr_anis[i] = 1;
	  cl11[i] = mxGetPr(prhs[4+num_fields+i+pofs]);
	  cl12[i] = mxGetPr(prhs[4+num_fields+i+pofs])+xsiz*ysiz;
	  cl22[i] = mxGetPr(prhs[4+num_fields+i+pofs])+2*xsiz*ysiz;
	} else
	  mexErrMsgTxt( "Wrong sized correlation length matrix." );
      }
      if (corr_spvar[i] > 0)
	if (!((corrtype[i] == 0) || (corrtype[i] == 2) || (corrtype[i] == 5)))
	  mexErrMsgTxt( "Nonstationary covariance for types 0,2,5 only." );
    } else {
      /* user has supplied explicit corr. fn */
      usercorrs[i] = mxGetPr(prhs[4+num_fields+i+pofs]);
      if ((mxGetM(prhs[4+num_fields+i+pofs]) % 2) + (mxGetN(prhs[4+num_fields+i+pofs]) % 2) < 2)
	mexErrMsgTxt( "Correlation matrix must be of size odd x odd." );
      userorigin[i][0] = (mxGetN(prhs[4+num_fields+i+pofs])-1)/2;
      userorigin[i][1] = (mxGetM(prhs[4+num_fields+i+pofs])-1)/2;
    }
    if (sigma_type[i])
      if ((mxGetN(prhs[4+i+pofs])!=xsiz) || (mxGetM(prhs[4+i+pofs])!=ysiz))
	mexErrMsgTxt( "Wrong sized standard deviation multipliers matrix." );
  }

  /* adjust correlation lengths per correl. function needs */
  for (f=0; f<num_fields; f++) {
    len_x[f] = correl_length( f, len_x[f] );
    len_y[f] = correl_length( f, len_y[f] );
  }

  /* field is stationary (unless space varying); code for local offsets */
  for (f=0; f<num_fields; f++)
    for (i=-LocalPix; i<LocalPix; i++)
      for (j=-LocalPix; j<LocalPix; j++) {
#ifdef Detailed_Distance
	localcorrs[f][LocalPix+i][LocalPix+j] = correl_stat(f,i,j,0);
#else
	localcorrs[f][LocalPix+i][LocalPix+j] = correl_stat(f,i,j);
#endif
	localdists[LocalPix+i][LocalPix+j] = sqrt(i*i+j*j);
      }

  return(4+2*num_fields+pofs);
}


/* Returns the modified Bessel function I1(x) for any real x. */
double bessi1(double x)
{
  double ax,ans;
  double y;
  if ((ax=fabs(x)) < 3.75) {
    y=x/3.75;
    y*=y;
    ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934+y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
  } else {
    y=3.75/ax;
    ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1-y*0.420059e-2));
    ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2+y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
    ans *= (exp(ax)/sqrt(ax));
  }
  return x < 0.0 ? -ans : ans;
}

/* Returns the modied Bessel function K1(x) for positive real x. */
double bessk1(double x)
{
  double bessi1(double x);
  double y,ans;
  if (x <= 2.0) {
    y=x*x/4.0;
    ans=(log(x/2.0)*bessi1(x))+(1.0/x)*(1.0+y*(0.15443144+y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1+y*(-0.110404e-2+y*(-0.4686e-4)))))));
  } else {
    y=2.0/x;
    ans=(exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619+y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2+y*(0.325614e-2+y*(-0.68245e-3)))))));
  }
  return ans;
}

int wrap_index( int f, int x, int y, int closest )
{
  /* needs to consider separate x,y wrapping */

  if (closest == 0) {
    /* don't want to find closest pixel, return good index or -1 */

    if (wraptype[f]) return( ysiz * ((x + xsiz) % xsiz) + ((y + ysiz) % ysiz));

    if ((x<0) || (y<0) || (x>=xsiz) || (y>=ysiz)) return(-1);
    return( ysiz * x + y );
  }
  else {
    /* if outside of region and not wrapping, find closest */

    if (wraptype[f]) return( ysiz * ((x + xsiz) % xsiz) + ((y + ysiz) % ysiz));
    return( max(min(y,ysiz-1),0) + ysiz * max(min(x,xsiz-1),0) );
  }
}

/* Main correlation function to be called externally */
#ifdef CorrFn_Region
double corrfn( int field1, int x1, int y1, int l1, int field2, int x2, int y2, int l2 )
{
  int i1, i2, d, ih;
  double v;
#else
double corrfn( int field1, int x1, int y1, int field2, int x2, int y2 )
{
  int i1, i2, d, ih;
  int l1=0, l2=0;
  double v;
#endif

  if (debug_mode && (debug_opt==2)) {

    printf( "Corr: F%d: (%3d,%3d)   F%d: (%3d,%3d)   ", field1,x1,y1,field2,x2,y2);
    printf( "Regions %1d,%1d  Index (%d,%d)\n", (int) landptr[y1+x1*ysiz],(int) landptr[y2+x2*ysiz], y1+x1*ysiz,y2+x2*ysiz );
  }

  if (field1 != field2) return( 0.0 );
  if (field1 >= num_fields) return( 0.0 );

  if ((x1==x2) && (y1==y2)) {
    if (!sigma_type[field1]) 
      return( sgm[field1][0] * sgm[field2][0] );
    
    return( sgm[field1][y1+x1*ysiz] * sgm[field2][y2+x2*ysiz] );
  }

  i1 = wrap_index(field1,x1,y1,0);
  i2 = wrap_index(field2,x2,y2,0);

  if (l1 == 0 && i1>=0) l1 = ((int) landptr[i1]);
  if (l2 == 0 && i2>=0) l2 = ((int) landptr[i2]);
  
  if ((l1 == 0) || (l2 == 0))
    return( 0.0 );

  if ((l1 <= adjsiz) && (l2 <= adjsiz)) {
    if (adjptr[field1][(l1-1)*adjsiz+l2-1]<0.5)
      return( 0.0 );
  }
  else if (l1 != l2)
    return( 0.0 );

  /* find closest index if field unwrapped and pixel outside */
  if (i1<0) i1 = wrap_index(field1,x1,y1,1);
  if (i2<0) i2 = wrap_index(field2,x2,y2,1);

#ifdef Detailed_Distance
  /* get local detailed distance */
  d = abs(distptr[i1]-distptr[i2]);
#endif

  if (!sigma_type[field1]) { i1 = i2 = 0; }

  /* make adjacent across boundaries left/right and top/bottom */
  if (wraptype[field1]) {
    l1 = ((x1-x2+xsiz+xsiz/2) % xsiz) - xsiz/2;
    l2 = ((y1-y2+ysiz+ysiz/2) % ysiz) - ysiz/2;
  }
  else {
    l1 = x1-x2; 
    l2 = y1-y2;
  }

  if (corr_spvar[field1] == 0) {
    if ((abs(l1) < LocalPix) && (abs(l2) < LocalPix))
#ifdef Detailed_Distance
      if (d <= localdists[LocalPix+l1][LocalPix+l2])
#endif
	return( sgm[field1][i1] * sgm[field2][i2] * localcorrs[field1][LocalPix+l1][LocalPix+l2] );

#ifdef Detailed_Distance
    return( sgm[field1][i1] * sgm[field2][i2] * correl_stat(field1,l1,l2,d) );
#else
    return( sgm[field1][i1] * sgm[field2][i2] * correl_stat(field1,l1,l2) );
#endif
  }

  /* find halfway point between pixels, wrap if necessary */
  ih = wrap_index( field1, x1+l1/2, y1+l2/2, 1 );

#ifdef Detailed_Distance
  return( sgm[field1][i1] * sgm[field2][i2] * correl_nonstat(field1,l1,l2,ih,d) );
#else
  return( sgm[field1][i1] * sgm[field2][i2] * correl_nonstat(field1,l1,l2,ih) );
#endif
}

double corrfn_quick( int field, int x1, int y1, int x2, int y2 )
{
  int i1, i2, d, ih;
  int l1=0, l2=0;
  double v;


  if (debug_mode && (debug_opt==2)) {
    printf( "QCorr: F%d: (%3d,%3d) (%3d,%3d)\n", field,x1,y1,x2,y2);
  }

  if ((x1==x2) && (y1==y2)) {
    if (!sigma_type[field]) {
      v = sgm[field][0]; return( v*v );
    }
    
    v = sgm[field][wrap_index(field,x1,y1,1)]; return( v*v );
  }

  /* make adjacent across boundaries left/right and top/bottom */
  if (wraptype[field]) {
    l1 = ((x1-x2+xsiz+xsiz/2) % xsiz) - xsiz/2;
    l2 = ((y1-y2+ysiz+ysiz/2) % ysiz) - ysiz/2;
  }
  else {
    l1 = x1-x2; 
    l2 = y1-y2;
  }

  i1 = i2 = 0;

#ifdef Detailed_Distance
  /* get local detailed distance */
  i1 = wrap_index(field,x1,y1,1);
  i2 = wrap_index(field,x2,y2,1);
  d = abs(distptr[i1]-distptr[i2]);
#else
  if (sigma_type[field]) {
    i1 = wrap_index(field,x1,y1,1);
    i2 = wrap_index(field,x2,y2,1);
  }
#endif

  if (debug_mode && (debug_opt==2)) {
    printf( "   (QCorr) F%d ind,l  (%3d,%3d) (%3d,%3d)\n", field,i1,l1,i2,l2);
  }

  if (corr_spvar[field] == 0) {
    if ((abs(l1) < LocalPix) && (abs(l2) < LocalPix))
#ifdef Detailed_Distance
      if (d <= localdists[LocalPix+l1][LocalPix+l2])
#endif
	return( sgm[field][i1] * sgm[field][i2] * localcorrs[field][LocalPix+l1][LocalPix+l2] );

#ifdef Detailed_Distance
    return( sgm[field][i1] * sgm[field][i2] * correl_stat(field,l1,l2,d) );
#else
    return( sgm[field][i1] * sgm[field][i2] * correl_stat(field,l1,l2) );
#endif
  }

  /* find halfway point between pixels, wrap if necessary */
  ih = wrap_index( field, x1+l1/2, y1+l2/2, 1 );

#ifdef Detailed_Distance
  return( sgm[field][i1] * sgm[field][i2] * correl_nonstat(field,l1,l2,ih,d) );
#else
  return( sgm[field][i1] * sgm[field][i2] * correl_nonstat(field,l1,l2,ih) );
#endif
}

double correl_length(int f, double l)
{
  /* want parameter as correlation length only for sep-exp */
  if (corrtype[f] == 1) 
    /* Separable Exponential */
    return(l);

  return(l*l);
}

double correl_fn( int f, double r )
{
  /* correction factors on r serve to normalize correlation lengths */
  /* e.g., factor of 1.88 in case 3, 1.66 in case 5 etc. */

  if (corrtype[f] == 0) {
    /* Gaussian */
    return( exp(-r) );
  }
  if (corrtype[f] == 2) {
    /* NonSeparable Exponential */
    return( exp(-sqrt(r)) );
  }
  if (corrtype[f] == 3) {
    /* Smooth Non-Gaussian */
    r = 1.88*r;
    return( (1-exp(-r))*exp(-sqrt(r))+exp(-r) );
  }
  if (corrtype[f] == 5) {
    /* Bessel */
    r = 1.66 * sqrt(r);
    if (r<0.0001) return(1.0);
    return( r*bessk1(r) );
  }
  if (corrtype[f] == 6) {
    /* Spherical */
    double rs;
    r = 0.45 * r;
    rs = sqrt(r);
    if (r<1) return( 1-1.5*rs+0.5*r*rs );
    return(0.0);
  }
  if (corrtype[f] == 7) {
    /* Logistic */
    r = r * 3.44;
    return( 1.0-r/(2.0+r) );
  }

  /* Default is Gaussian */
  return( exp(-r) );
}


#ifdef Detailed_Distance
double correl_stat( int f, int i, int j, double d )
#else
double correl_stat( int f, int i, int j )
#endif
{
  double ds, ijs;

  /* two special cases to handle here, the rest are standard */
  if (corrtype[f] == 1) {
    /* Separable Exponential */
#ifdef Detailed_Distance
    if (d*d > i*i+j*j) { 
      double dr; 
      dr = d/sqrt(i*i+j*j);
      return( exp(-dr*abs(i)/len_x[f]-dr*abs(j)/len_y[f]) );
    }
#endif
    return( exp(-abs(i)/len_x[f]-abs(j)/len_y[f]) );
  }
  if (corrtype[f] == 4) {
    /* User-supplied correlation */
#ifdef Detailed_Distance
    if (d*d > i*i+j*j) { 
      double dr; 
      dr = d/sqrt(i*i+j*j);
      i = (int) round(dr*i);
      j = (int) round(dr*j);
    }
#endif
    if ((abs(i)<=userorigin[f][0]) && (abs(j)<=userorigin[f][1]))
      return( usercorrs[f][(userorigin[f][0]+i)*(userorigin[f][1]*2+1)+userorigin[f][1]+j] );
    return( 0.0 );
  }

  /* continue with standard case */
#ifdef Detailed_Distance
  ds = d*d; ijs=i*i+j*j;
  if (ds > ijs) 
    return( correl_fn( f, ds*i*i/(len_x[f]*ijs) + ds*j*j/(len_y[f]*ijs) ) );
#endif
  return( correl_fn( f, i*i/len_x[f] + j*j/len_y[f] ) );
}

//static nswarning=0;
static int nswarning=0;
#ifdef Detailed_Distance
double correl_nonstat( int f, int i, int j, int k, double d )
#else
double correl_nonstat( int f, int i, int j, int k )
#endif
{
  int is,js,ds;
  double cl;

  if (nswarning==0) {
    printf("Warning -- Nonstationary statistics selected.\n");
    nswarning=1;
  }

  is = i*i; js=j*j;
#ifdef Detailed_Distance
  ds = d*d;
  if (ds > is+js) { 
    if (corr_anis[f]==0) {
      cl = cl11[f][k];
      return( correl_fn( f, ds/(cl*cl) ) );
    } else {
      return( correl_fn( f, ds*(is*cl11[f][k]+2*i*j*cl12[f][k]+js*cl22[f][k])/(is+js) ) );
    }
  }
#endif
  if (corr_anis[f]==0) {
    cl = cl11[f][k];
    return( correl_fn( f, (i*i+j*j)/(cl*cl) ) );
  } else {
    return( correl_fn( f, i*i*cl11[f][k]+2*i*j*cl12[f][k]+j*j*cl22[f][k] ) );
  }
}



/* helper routines to determine pixel regions */
int corr_region( int f, int x, int y, int wx, int wy )
{
  /* return:  0    - all null, 
              num  - all one group, possibly plus null
              -num - multiple groups */

  int grp, ind, dx, dy, l;

  grp = 0;
  for (dx=x;dx<x+wx;dx++) 
    for (dy=y;dy<y+wy;dy++) {
      ind = wrap_index( f,dx,dy,0 );
      if (ind >= 0) {
	if ((l = ((int) landptr[ind])) > 0) {
	  if (grp == 0) grp = l; 
	  if (l != grp) {
	    if (l <= adjsiz) {
	      if (adjptr[f][(l-1)*adjsiz+grp-1]<0.5) return( -1 );
	    } else
	      return( -1 );
	  }
	}
      }
    }
  return( grp );
}

int corr_pair( int f, int x1, int y1, int x2, int y2 )
{
  int i1, i2, l1, l2;

  if ((x1==x2) && (y1==y2)) return( 1 );

  i1 = wrap_index( f,x1,y1,0 );
  i2 = wrap_index( f,x2,y2,0 );
  if (i1<0 || i2<0) return(0);

  l1 = ((int) landptr[i1]);
  l2 = ((int) landptr[i2]);
  
  if ((l1 == 0) || (l2 == 0)) return( 0 );

  if ((l1 <= adjsiz) && (l2 <= adjsiz)) {
    if (adjptr[f][(l1-1)*adjsiz+l2-1]<0.5)
      return( 0 );
  }
  else if (l1 != l2)
    return( 0 );

  return( 1 );
}

int corr_types( int f, int l1, int l2 )
{
  if ((l1 == 0) || (l2 == 0)) return( 0 );

  if ((l1 <= adjsiz) && (l2 <= adjsiz)) {
    if (adjptr[f][(l1-1)*adjsiz+l2-1]<0.5)
      return( 0 );
  }
  else if (l1 != l2)
    return( 0 );

  return( 1 );
}

int corr_null( int f, int x, int y )
{
  int i;

  if ((i = wrap_index(f,x,y,0))<0) return(1);

  if (((int) landptr[i]) == 0) return(1);
  return(0);
}

int corr_type( int f, int x, int y )
{
  int i;

  if ((i = wrap_index(f,x,y,0))<0) return(0);

  return( (int) landptr[i] );
}

