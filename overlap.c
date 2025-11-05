/**************************************************************************
  
   overlap.c

   This function performs a variety of utility operation on 
   matrices for the purpose of working with trees containing
   overlapping representations.  The function is expected
   to be used from within the Matlab environment.

   Compiler variables that may be specified by the user:
      Print_Comments    - if defined, enables printing of comments


   Compilation:   mex overlap.c


   out_mat = overlap( in_mat, num_scales, dx, dy, function )
     in_mat     - a dense or overlapped matrix, depending on function
     overlap    - region overlap as a function of scale (fine ... coarse)
     dx         - tree splitting from one scale to next along x
     dy         - tree splitting from one scale to next along y
     function   - operation to perform:
                    0 - copy values from in_mat to overlapped out_mat
                    1 - sum values from overlapped in_mat to dense out_mat
                    2 - determine linear weighting coefficients for overlap
                        and write to out_mat (obsolete)
                    3 - weights, using one-sided sloping (recommended)

   Copyright (c) Paul Fieguth
   April 1994

   Dec, 1999  - Allow dx,dy to be vectors, not just scalars

*/


#include <math.h>
#include "mex.h"

#define Max_Scales 12

#define	max(A, B)	((A) > (B) ? (A) : (B))
#define	min(A, B)	((A) < (B) ? (A) : (B))


#define M_mat      prhs[0]
#define M_ovlp     prhs[1]
#define M_dx       prhs[2]
#define M_dy       prhs[3]
#define M_function prhs[4]
#define function ((int) floor(mxGetScalar( M_function )+0.5))


/* Global variables (required by recursive routines below main gateway) */
double *p_mat, *p_outmat, *p_ovlp, *p_weights[Max_Scales];
int osizex, osizey, sizex, sizey;
int dx[Max_Scales], dy[Max_Scales], scales;

/* Recursive routine function prototypes */
void copy_to_overlap( int scale, int lox, int hix, int loy, int hiy, 
		     int ix, int iy, int wx, int wy );
void sum_to_dense( int scale, int lox, int hix, int loy, int hiy, 
		  int ix, int iy, int wx, int wy );
void determine_overlap_weights( int scale, int lox, int hix, int loy, int hiy, 
			       int ix, int iy, int wx, int wy );
void determine_onesided_weights( int scale, int lox, int hix, int loy, int hiy,
			       int ix, int iy, int wx, int wy );


/*
 * Gateway routine
 */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  double temp_double;
  int i, j, k;
  
  /* Check for proper number of arguments */
  if (nrhs != 5) {
    mexPrintf( "Correct function usage:\n" );
    mexPrintf( "overlap( in_mat, overlap, dx, dy, function )\n" );
    mexPrintf( "   in_mat     - a dense or overlapped matrix, depending on function\n" );
    mexPrintf( "   overlap  - region overlap as a function of scale (coarse ... fine)\n" );
    mexPrintf( "   dx       - tree splitting from one scale to next along x\n" );
    mexPrintf( "   dy       - tree splitting from one scale to next along y\n" );
    mexPrintf( "   function - operation to perform:\n" );
    mexPrintf( "        0 - copy values from in_mat to overlapped out_mat\n" );
    mexPrintf( "        1 - sum values from overlapped in_mat to dense out_mat\n" );
    mexPrintf( "        3 - determine linear weighting coefficients for overlap\n" );
    mexPrintf( "            and write to out_mat (in_mat is ignored)\n" );
    mexErrMsgTxt("OVERLAP requires five input arguments.\n");
  }

  /* copy globally required parameters */
  scales = mxGetM(M_ovlp)+1;

  /* Check the dimensions of inputs */
  if ((scales < 1) || (scales > Max_Scales))
    mexErrMsgTxt( "Invalid number of scales implied by overlap vector" );
  if (mxGetN(M_ovlp) != 2)
    mexErrMsgTxt( "Overlap matrix must have two columns" );
  if ((function < 0) || (function > 3))
    mexErrMsgTxt( "Invalid function option specified" );
  i = mxGetM(M_dx)*mxGetN(M_dx);
  j = mxGetM(M_dy)*mxGetN(M_dy);
  if ((i != 1) && (i != (scales-1))) mexErrMsgTxt( "Invalid parameter dx" );
  if ((j != 1) && (j != (scales-1))) mexErrMsgTxt( "Invalid parameter dx" );
  for (k=0; k<scales-1; k++) {
    dx[k] = (int) floor(mxGetPr(M_dx)[(i==1)?0:k]+0.5);
    dy[k] = (int) floor(mxGetPr(M_dy)[(j==1)?0:k]+0.5);
  }

  /* link pointers to matrices */
  p_mat = mxGetPr( M_mat );
  p_ovlp = mxGetPr( M_ovlp );

  /* infer matrix sizes */
  osizex = osizey = 1;
  for (i=1; i<scales; i++) { osizex *= dx[i-1]; osizey *= dy[i-1]; }

  sizex = sizey = 1;
  for (i=scales-2; i>=0; i--) { 
    sizex = sizex * dx[i] - (dx[i]-1) * p_ovlp[i];
    sizey = sizey * dy[i] - (dy[i]-1) * p_ovlp[scales-1+i];
  }

#ifdef Print_Comments
  printf( "Overlap:  Scales %d, Ovx %d, Ovy %d, Denx %d, Deny %d\n", 
	   scales, osizex, osizey, sizex, sizey );
#endif

  if (function == 0) {
    /*********************************************************************
     * 
     * copy values from dense form into multiple overlapped regions 
     *
     */

    /* test input matrix */
    if ((mxGetM(M_mat) != sizey) || (mxGetN(M_mat) != sizex))
      mexErrMsgTxt( "Invalid dense matrix dimensions." );

    /* create overlapped output matrix */
    plhs[0] = mxCreateDoubleMatrix( osizey, osizex, mxREAL );
    p_outmat = mxGetPr( plhs[0] );

    /* call recursive copy routine */
    copy_to_overlap( 0, 0, sizex-1, 0, sizey-1, 0, 0, osizex, osizey );
  }
  else if (function == 1) {
    /*********************************************************************
     * 
     * sum values from overlapped matrix into dense form 
     *
     */

    /* test input matrix */
    if ((mxGetM(M_mat) != osizey) || (mxGetN(M_mat) != osizex))
      mexErrMsgTxt( "Invalid overlap matrix dimensions" );

    /* create dense output matrix */
    plhs[0] = mxCreateDoubleMatrix( sizey, sizex, mxREAL );
    p_outmat = mxGetPr( plhs[0] );

    /* clear matrix */
    for (i=0; i<sizex; i++) for (j=0; j<sizey; j++)
      p_outmat[j+sizey*i] = 0.0;

    /* call recursive copy routine */
    sum_to_dense( 0, 0, sizex-1, 0, sizey-1, 0, 0, osizex, osizey );
  }
  else if (function >= 2) {
    /*********************************************************************
     * 
     * determine weights for merging overlapped values into dense form
     *
     */

    /* create overlapped output matrix */
    plhs[0] = mxCreateDoubleMatrix( osizey, osizex, mxREAL );
    p_outmat = mxGetPr( plhs[0] );
    for (i=0; i<osizex; i++) for (j=0; j<osizey; j++)
      p_outmat[j+osizey*i] = -1;

    /* create weight matrices for each scale */
    sizex = sizey = 1;
    p_weights[scales-1] = &temp_double;
    for (i=scales-2; i>=0; i--) { 
      sizex = sizex * dx[i] - (dx[i]-1) * p_ovlp[i];
      sizey = sizey * dy[i] - (dy[i]-1) * p_ovlp[scales-1+i];
      p_weights[i] = mxCalloc( sizex * sizey, sizeof( double ) );
    }

    /* initialize weights at coarsest level */
    for (i=0; i<sizex; i++) for (j=0; j<sizey; j++)
      p_weights[0][j+i*sizey] = 1.0;
    
    /* call recursive copy routine */
    if (function == 2) 
      determine_overlap_weights( 0, 0, sizex-1, 0, sizey-1, 0, 0, osizex, osizey );
    if (function == 3)
      determine_onesided_weights( 0, 0, sizex-1, 0, sizey-1, 0, 0, osizex, osizey );
    /* free memory */
    for (i=scales-2; i>=0; i--)
      mxFree( p_weights[i] );
  }  
  return;
}


void copy_to_overlap( int scale, int lox, int hix, int loy, int hiy, 
		     int ix, int iy, int wx, int wy )
{
  int x, y, sublox, subloy, lenx, leny;

  if (scale == scales-1) {
    /* at finest level -- just copy appropriate value */
    p_outmat[iy+ix*osizey] = p_mat[loy+sizey*lox];
  }
  else {
    /* not at finest level, split and descend */
    for (x=0; x<dx[scale]; x++) for (y=0; y<dy[scale]; y++) {
      lenx = (hix-lox+1+(dx[scale]-1)*p_ovlp[scale])/dx[scale];
      leny = (hiy-loy+1+(dy[scale]-1)*p_ovlp[scales-1+scale])/dy[scale];
      sublox = lox+x*(lenx-p_ovlp[scale]);
      subloy = loy+y*(leny-p_ovlp[scales-1+scale]);
      copy_to_overlap( scale+1, sublox, sublox+lenx-1, subloy, subloy+leny-1, 
		      ix + x*wx/dx[scale], iy + y*wy/dy[scale], 
		      wx/dx[scale], wy/dy[scale] );
    }
  }
  return;
}

void sum_to_dense( int scale, int lox, int hix, int loy, int hiy, 
		   int ix, int iy, int wx, int wy )
{
  int x, y, sublox, subloy, lenx, leny;

  if (scale == scales-1) {
    /* at finest level -- just add appropriate value */
    p_outmat[loy+sizey*lox] += p_mat[iy+ix*osizey];
  }
  else {
    /* not at finest level, split and descend */
    for (x=0; x<dx[scale]; x++) for (y=0; y<dy[scale]; y++) {
      lenx = (hix-lox+1+(dx[scale]-1)*p_ovlp[scale])/dx[scale];
      leny = (hiy-loy+1+(dy[scale]-1)*p_ovlp[scales-1+scale])/dy[scale];
      sublox = lox+x*(lenx-p_ovlp[scale]);
      subloy = loy+y*(leny-p_ovlp[scales-1+scale]);
      sum_to_dense( scale+1, sublox, sublox+lenx-1, subloy, subloy+leny-1, 
		   ix + x*wx/dx[scale], iy + y*wy/dy[scale], 
		   wx/dx[scale], wy/dy[scale] );
    }
  }
  return;
}

void determine_overlap_weights( int scale, int lox, int hix, int loy, int hiy, 
			       int ix, int iy, int wx, int wy )
{
  int x, y, sublox, subloy, lenx, leny;
  int i, j;

  if (scale == scales-1) {
    /* at finest level, just copy weight to output matrix */
    p_outmat[iy+ix*osizey] = p_weights[scale][0];
  }
  else {
    /* not at finest level - split, set weights, and descend */
    /* not at finest level, split and descend */
    for (x=0; x<dx[scale]; x++) for (y=0; y<dy[scale]; y++) {
      lenx = (hix-lox+1+(dx[scale]-1)*p_ovlp[scale])/dx[scale];
      leny = (hiy-loy+1+(dy[scale]-1)*p_ovlp[scales-1+scale])/dy[scale];

      sublox = x*(lenx-p_ovlp[scale]);
      subloy = y*(leny-p_ovlp[scales-1+scale]);

      /* set weights in this child matrix */
      for (i=0; i<lenx; i++) for (j=0; j<leny; j++)
	p_weights[scale+1][j+leny*i] = 
	  p_weights[scale][(subloy+j)+(sublox+i)*(hiy-loy+1)];

      for (i=0; i<lenx; i++) for (j=0; j<p_ovlp[scales-1+scale]; j++) {
	p_weights[scale+1][j+leny*i] *= (j+1)/(p_ovlp[scales-1+scale]+1);
	p_weights[scale+1][leny-j-1+leny*i] *= (j+1)/(p_ovlp[scales-1+scale]+1);
      }

      for (i=0; i<leny; i++) for (j=0; j<p_ovlp[scale]; j++) {
	p_weights[scale+1][i+leny*j] *= (j+1)/(p_ovlp[scale]+1);
	p_weights[scale+1][i+leny*(lenx-1-j)] *= (j+1)/(p_ovlp[scale]+1);
      }
	
      /* descend to next scale */
      sublox += lox;
      subloy += loy;
      determine_overlap_weights( scale+1, sublox, sublox+lenx-1, subloy, subloy+leny-1, 
				ix + x*wx/dx[scale], iy + y*wy/dy[scale], 
				wx/dx[scale], wy/dy[scale] );
    }
  }    
}

void determine_onesided_weights( int scale, int lox, int hix, int loy, int hiy,
			       int ix, int iy, int wx, int wy )
{
  int x, y, sublox, subloy, lenx, leny;
  int i, j;

  if (scale == scales-1) {
    /* at finest level, just copy weight to output matrix */
    p_outmat[iy+ix*osizey] = p_weights[scale][0];
  }
  else {
    /* not at finest level - split, set weights, and descend */
    for (x=0; x<dx[scale]; x++) for (y=0; y<dy[scale]; y++) {
      lenx = (hix-lox+1+(dx[scale]-1)*p_ovlp[scale])/dx[scale];
      leny = (hiy-loy+1+(dy[scale]-1)*p_ovlp[scales-1+scale])/dy[scale];

      sublox = x*(lenx-p_ovlp[scale]);
      subloy = y*(leny-p_ovlp[scales-1+scale]);

      /* set weights in this child matrix */
      for (i=0; i<lenx; i++) for (j=0; j<leny; j++)
	p_weights[scale+1][j+leny*i] = 
	  p_weights[scale][(subloy+j)+(sublox+i)*(hiy-loy+1)];

      for (i=0; i<lenx; i++) for (j=0; j<p_ovlp[scales-1+scale]; j++) {
	if (y > 0)
	  p_weights[scale+1][j+leny*i] *= (j+1)/(p_ovlp[scales-1+scale]+1);
	if (y < dy[scale]-1)
	  p_weights[scale+1][leny-j-1+leny*i] *= (j+1)/(p_ovlp[scales-1+scale]+1);
      }

      for (i=0; i<leny; i++) for (j=0; j<p_ovlp[scale]; j++) {
	if (x > 0)
	  p_weights[scale+1][i+leny*j] *= (j+1)/(p_ovlp[scale]+1);
	if (x < dx[scale]-1)
	  p_weights[scale+1][i+leny*(lenx-1-j)] *= (j+1)/(p_ovlp[scale]+1);
      }
	
      /* descend to next scale */
      sublox += lox;
      subloy += loy;
      determine_onesided_weights( scale+1, sublox, sublox+lenx-1, subloy, subloy+leny-1, 
				 ix + x*wx/dx[scale], iy + y*wy/dy[scale], 
				 wx/dx[scale], wy/dy[scale] );
    }
  }    
}

