#include <stdio.h>
#include "smth.h"

/* this file computes the effective position of the passed node in space.
   This is useful for moving anisotropies from physical space onto a
   multiscale model in the overlapped domain.  The function returns the
   location of the corner of the domain, and the size of the domain in
   finest-scale pixel units */

#define MaxScales 20

int ovlpinit = 0;
int *ovlpmem = NULL;
int *ovlpxposn[MaxScales], *ovlpyposn[MaxScales];
int *ovlpxsiz[MaxScales], *ovlpysiz[MaxScales];
int ovlpynodes[MaxScales];

void overlap_init( int scales, int *ovlp, int *dx, int *dy )
{
  int scale, i, t1, t2, sx, sy;

  ovlpinit = 0;
  sx = 1; sy = 1;
  pf_Malloc( (void **) &ovlpmem, sizeof(int), 0 );

  for (scale=0; scale<scales; scale++) {
    pf_Malloc( (void **) ovlpxposn+scale, sx*sizeof(int), 0 );
    pf_Malloc( (void **) ovlpyposn+scale, sy*sizeof(int), 0 );
    pf_Malloc( (void **) ovlpxsiz+scale, sx*sizeof(int), 0 );
    pf_Malloc( (void **) ovlpysiz+scale, sy*sizeof(int), 0 );
    ovlpynodes[scale] = sy;

    for (i=0; i<sy; i++) 
      overlap_data( scale, i, scales, ovlp, dx, dy, &t1, ovlpyposn[scale]+i, &t2, ovlpysiz[scale]+i );

    for (i=0; i<sx; i++) 
      overlap_data( scale, sy*i, scales, ovlp, dx, dy, ovlpxposn[scale]+i, &t1, ovlpxsiz[scale]+i, &t2 );

    sx *= dx[scale]; sy *= dy[scale];
  }
  ovlpinit = 1;
}

void overlap_data( int scale, int s, int scales, int *ovlp, int *dx, int *dy, int *x, int *y, int *xsiz, int *ysiz )
{
  int sxsiz, sysiz, sx, sy, ix, iy, i;

  if (ovlpmem == NULL)
    overlap_init( scales, ovlp, dx, dy );

  if (ovlpinit) {
    ix = s/ovlpynodes[scale];
    iy = s%ovlpynodes[scale];
    *x = ovlpxposn[scale][ix];
    *y = ovlpyposn[scale][iy];
    *xsiz = ovlpxsiz[scale][ix];
    *ysiz = ovlpysiz[scale][iy];
    return;
  }

  /* initialize position variables */
  *x = *y = 0;

  /* initialize sizing variables */
  *xsiz = 1;
  *ysiz = 1;

  for (i=scales-2; i>=0; i--) {
    *xsiz = *xsiz * dx[i] - (dx[i]-1)*ovlp[i];
    *ysiz = *ysiz * dy[i] - (dy[i]-1)*ovlp[i+scales-1];
  }

  sxsiz = sysiz = 1;
  for (i=0; i<scale; i++) {
    sxsiz *= dx[i];
    sysiz *= dy[i];
  }
  sx = s / sysiz;
  sy = s-sx*sysiz;

  /* loop down to scale and infer position */
  for (i=0; i<scale; i++) {
    ix = sx/(sxsiz/dx[i]);
    iy = sy/(sysiz/dy[i]);
    
    *xsiz = (*xsiz+(dx[i]-1)*ovlp[i])/dx[i];
    *ysiz = (*ysiz+(dy[i]-1)*ovlp[i+scales-1])/dy[i];

    *x += ix * (*xsiz-ovlp[i]);
    *y += iy * (*ysiz-ovlp[i+scales-1]);

    sxsiz /= dx[i];
    sysiz /= dy[i];

    sx -= ix * sxsiz;
    sy -= iy * sysiz;
  }

}
