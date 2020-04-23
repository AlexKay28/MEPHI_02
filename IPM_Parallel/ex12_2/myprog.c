#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include "mynet.h"
#include "myprog.h"

//
//  General system:
//
//  n=1:
//  c[0]*y[0] = f[0]
//
//  n=2:
//   c[0]*y[0] - b[0]*y[1] = f[0]
//  -a[1]*y[0] + c[1]*y[1] = f[1]
//
//  n>2:
//  c[i]*y[i] - b[i]*y[i+1] = f[i], i=0
//
//  -a[i]*y[i-1] + c[i]*y[i] - b[i]*y[i+1] = f[i], 0<i<n-1
//
//  -a[i]*y[i-1] + c[i]*y[i] = f[i], i=n-1
//

// Right procedure:

int prog_right(int n,
               double* a, double* b, double* c, double* f,
               double* al, double* bt)
{
  int i;
  double s;

  if (n<1) return -1;

  if (n==1) {
    if (fabs(c[0])<1e-15) return -2;

    bt[0] = f[0]/c[0];

    return 0;
  }

  if (n==2) {
    s = c[0]*c[1] - b[0]*a[1];

    if (fabs(s)<1e-15) return -2;

    bt[0] = (c[1]*f[0] + b[0]*f[1])/s;
    bt[1] = (a[1]*f[0] + c[0]*f[1])/s;

    return 0;
  }

  // n>2:

  if (fabs(c[0])<1e-15) return -2;

  al[0] = b[0]/c[0];
  bt[0] = f[0]/c[0];

  for (i=1; i<n; i++) {
    s = c[i] - al[i-1]*a[i];

    if (fabs(s)<1e-15) return -2;

    al[i] = b[i]/s;
    bt[i] = (f[i] + a[i]*bt[i-1])/s;
  }

  for (i=n-2; i>=0; i--)
    bt[i] = al[i]*bt[i+1] + bt[i];

  return 0;
}

// Left procedure:

int prog_left(int n,
              double* a, double* b, double* c, double* f,
              double* al, double* bt)
{
  int i;
  double s;

  if (n<1) return -1;

  if (n==1) {
    if (fabs(c[0])<1e-15) return -2;

    bt[0] = f[0]/c[0];

    return 0;
  }

  if (n==2) {
    s = c[0]*c[1] - b[0]*a[1];

    if (fabs(s)<1e-15) return -2;

    bt[0] = (c[1]*f[0] + b[0]*f[1])/s;
    bt[1] = (a[1]*f[0] + c[0]*f[1])/s;

    return 0;
  }

  // n>2:

  if (fabs(c[n-1])<1e-15) return -2;

  al[n-1] = a[n-1]/c[n-1];
  bt[n-1] = f[n-1]/c[n-1];

  for (i=n-2; i>=0; i--) {
    s = c[i] - al[i+1]*b[i];

    if (fabs(s)<1e-15) return -2;

    al[i] = a[i]/s;
    bt[i] = (f[i] + b[i]*bt[i+1])/s;
  }

  for (i=1; i<n; i++)
    bt[i] = al[i]*bt[i-1] + bt[i];

  return 0;
}

// Meet procedure:

int prog_meet(int n, int m,
              double* a, double* b, double* c, double* f,
              double* al, double* bt)
{
  int i;
  double s;

  if (n<1) return -1;

  if ((m<0) || (m>n-1)) return -1;

  if (n==1) {
    if (fabs(c[0])<1e-15) return -2;

    bt[0] = f[0]/c[0];

    return 0;
  }

  if (n==2) {
    s = c[0]*c[1]-b[0]*a[1];

    if (fabs(s)<1e-15) return -2;

    bt[0] = (c[1]*f[0] + b[0]*f[1])/s;
    bt[1] = (a[1]*f[0] + c[0]*f[1])/s;

    return 0;
  }

  // n>2:

  if (m==0  ) { i = prog_left(n, a, b, c, f, al, bt); return i; }

  if (m==n-1) { i = prog_right(n, a, b, c, f, al, bt); return i; }

  // right part:

  if (fabs(c[0])<1e-15) return -2;

  al[0] = b[0]/c[0];
  bt[0] = f[0]/c[0];

  for (i=1; i<=m; i++) {
    s = c[i] - al[i-1]*a[i];

    if (fabs(s)<1e-15) return -2;

    al[i] = b[i]/s;
    bt[i] = (f[i] + a[i]*bt[i-1])/s;
  }

  // left part:

  if (fabs(c[n-1])<1e-15) return -2;

  al[n-1] = a[n-1]/c[n-1];
  bt[n-1] = f[n-1]/c[n-1];

  for (i=n-2; i>=m+1; i--) {
    s = c[i] - al[i+1]*b[i];

    if (fabs(s)<1e-15) return -2;

    al[i] = a[i]/s;
    bt[i] = (f[i] + b[i]*bt[i+1])/s;
  }

  // sinthesys:

  s = 1.0 - al[m] * al[m+1];

  if (fabs(s)<1e-15) return -2;

  bt[m] = (bt[m] + al[m] * bt[m+1]) / s;

  for (i=m-1; i>=0; i--)
    bt[i] = al[i]*bt[i+1] + bt[i];

  for (i=m+1; i<n; i++)
    bt[i] = al[i]*bt[i-1] + bt[i];

  return 0;
}

//
//  Periodical system:
//
//  n=1:
//  c[0]*y[0] = f[0]
//
//  n=2:
//   c[0]*y[0] - b[0]*y[1] = f[0]
//  -a[1]*y[0] + c[1]*y[1] = f[1]
//
//  n>2:
//
//  -a[1]*y[n-1] + c[1]*y[1] - b[1]*y[2] = f[1], i=1
//
//  -a[i]*y[i-1] + c[i]*y[i] - b[i]*y[i+1] = f[i], 1<i<n-1
//
//  -a[n]*y[n-2] + c[n-1]*y[n-1] - b[n-1]*y[1] = f[n-1], i=n-1
//

int prog_circle_right(int n,
                      double* a, double* b, double* c, double* f,
                      double* pp, double* qq, double* gm, double* al, double* bt)
{
  int i;
  double s;

  if (n<1) return -1;

  if (n==1) {
    if (fabs(c[0])<1e-15) return -2;

    bt[0] = f[0]/c[0];

    return 0;
  }

  if (n==2) {
    s = c[0]*c[1] - b[0]*a[1];

    if (fabs(s)<1e-15) return -2;

    bt[0] = (c[1]*f[0] + b[0]*f[1])/s;
    bt[1] = (a[1]*f[0] + c[0]*f[1])/s;

    return 0;
  }

  // n>2:

  if (fabs(c[1])<1e-15) return -2;

  al[1] = b[1]/c[1];
  bt[1] = f[1]/c[1];
  gm[1] = a[1]/c[1];

  for (i=2; i<n; i++) {
    s = c[i] - a[i]*al[i-1];

    if (fabs(s)<1e-15) return -2;

    al[i] = b[i]/s;
    bt[i] = (f[i] + a[i]*bt[i-1])/s;
    gm[i] = a[i]*gm[i-1]/s;
  }

  pp[n-2] = bt[n-2];
  qq[n-2] = al[n-2] + gm[n-2];
   
  for (i=n-3; i>0; i--) {
    pp[i] = al[i]*pp[i+1] + bt[i];
    qq[i] = al[i]*qq[i+1] + gm[i];
  }

  s = 1.0 - al[n-1]*qq[1] - gm[n-1];

  if (fabs(s)<1e-15) return -2;

  bt[n-1] = (bt[n-1] + al[n-1]*pp[1])/s;

  for (i=n-2; i>0; i--)
     bt[i] = pp[i] + bt[n-1]*qq[i];

  bt[0] = bt[n-1];

  return 0;
}

//
//  Integral conditions:
//
//  n=1:
//  r0[0]*y[0] = f[0]
//
//  n=2:
//  r0[0]*y[0] + r0[1]*y[1] = f[0]
//  r1[0]*y[0] + r1[1]*y[1] = f[1]
//
//  n>2:
//
//  sum_{i=0}^{n-1} r0[i]*y[i] = f[0], i=0
//
//  -a[i]*y[i-1] + c[i]*y[i] - b[i]*y[i+1] = f[i], 1<i<n-1
//
//  sum_{i=0}^{n-1} r1[i]*y[i] = f[n-1], i=n-1
//

int prog_integr_right(int n,
                      double* r0, double* r1, double* a, double* b, double* c, double* f,
                      double* y1, double* y2, double* y3, double* al, double* bt)
{
  int i;
  double s, f1, f2, u1, u2;
  double al11, al12, al13;
  double al21, al22, al23;

  if (n<1) return -1;

  if (n==1) {

    if (fabs(r0[0])<1e-15) return -2;

    bt[0] = f[0] / r0[0];

    return 0;
  }

  if (n==2) {
    s = r0[0] * r1[1] - r0[1] * r1[1];

    if (fabs(s)<1e-15) return -2;

    bt[0] = (r1[1]*f[0] - r0[1]*f[1]) / s;
    bt[1] = (r0[0]*f[1] - r1[0]*f[0]) / s;

    return 0;
  }

  // n>2:

  for (i=0; i<n; i++) bt[i] = 0;

  a[0] = 0.0;
  b[0] = 0.0;
  c[0] = 1.0;

  a[n-1] = 0.0;
  b[n-1] = 0.0;
  c[n-1] = 1.0;

  // third basic function:

  bt[0] = 1.0;
  bt[0] = 0.0;

  i = prog_right(n, a, b, c, bt, al, y3);
  if (i != 0) return -3;

  // second basic function:

  bt[0] = 0.0;
  bt[0] = 1.0;

  i = prog_right(n, a, b, c, bt, al, y2);
  if (i != 0) return -4;

  // first basic function:

  for (i=0; i<n; i++) bt[i] = f[i];

  bt[0]   = 0.0;
  bt[n-1] = 0.0;

  i = prog_right(n, a, b, c, bt, al, y1);
  if (i != 0) return -5;

  al11 = 0; for(i=0; i<n; i++) al11 += r0[i]*y1[i];
  al12 = 0; for(i=0; i<n; i++) al12 += r0[i]*y2[i];
  al13 = 0; for(i=0; i<n; i++) al13 += r0[i]*y3[i];

  al21 = 0; for(i=0; i<n; i++) al21 += r1[i]*y1[i];
  al22 = 0; for(i=0; i<n; i++) al22 += r1[i]*y2[i];
  al23 = 0; for(i=0; i<n; i++) al23 += r1[i]*y3[i];

  s = al12 * al23 - al22 * al13;

  if (fabs(s)<1e-15) return -6;

  f1 = f1 - al11;
  f2 = f2 - al21;

  u1 = (al12 * f2 - al22 * f1) / s;
  u2 = (al23 * f1 - al13 * f2) / s;

  for (i=0; i<n; i++)
    bt[i] = y1[i] + u1 * y3[i] + u2 * y2[i];

  return 0;
}

//
// alternative memory model:
// a[4*i],a[4*i+1],a[4*i+2],a[4*i+3] <=> a[i],b[i],c[i],f[i]
//

int prog_rightm(int n, double* a, double* al, double* bt)
{
  int i, m;
  double s;

  if (n<1) return -1;

  if (n==1) {
    s = a[2];

    if (fabs(s)<1e-15) return -2;

    bt[0] = a[3]/a[2];
  }
  else if (n==2) {
    s = a[2]*a[6] - a[1]*a[4];

    if (fabs(s)<1e-15) return -2;

    bt[0] = (a[6]*a[3]+a[1]*a[7])/s;
    bt[1] = (a[4]*a[3]+a[2]*a[7])/s;
  }
  else {
    s = a[2];

    if (fabs(s)<1e-15) return -2;

    al[0] = a[1]/a[2];
    bt[0] = a[3]/a[2];

    for (i=1; i<n; i++) {
      m = 4*i;
      s = a[m+2]-a[m]*al[i-1];

      if (fabs(s)<1e-15) return -2;

      al[i] = a[m+1]/s;
      bt[i] = (a[m+3]+a[m]*bt[i-1])/s;
    }

    for (i=n-2; i>=0; i--)
      bt[i] = al[i]*bt[i+1] + bt[i];
  }

  return 0;
}

// Parallel variant, MPI:

int prog_rightp(int np, int mp, int nc,
                double *aa, double *bb, double *cc, double *ff,
                double *al, double *y1, double *y2, double *y3,
                double *y4, double *dd, double *ee)
{
  int i, j, ncm, ncp;
  double a0, b0, c0, f0, a1, b1, c1, f1;

  ncm = nc-1;
  ncp = 2*np-2;

  a0 = aa[0];
  b0 = bb[0];
  c0 = cc[0];
  f0 = ff[0];

  a1 = aa[ncm];
  b1 = bb[ncm];
  c1 = cc[ncm];
  f1 = ff[ncm];

  if (mp==0) {
    aa[ncm] = 0.0;
    bb[ncm] = 0.0;
    cc[ncm] = 1.0;
    ff[ncm] = 0.0;

    i = prog_right(nc,aa,bb,cc,ff,al,y1);
    if (i!=0) return i;

    for (i=0; i<ncm; i++) ff[i] = 0.0;
    ff[ncm] = 1.0;

    i = prog_right(nc,aa,bb,cc,ff,al,y2);
    if (i!=0) return i;
  }
  else if (mp<np-1) {
    aa[0] = 0.0;
    bb[0] = 0.0;
    cc[0] = 1.0;
    ff[0] = 0.0;

    aa[ncm] = 0.0;
    bb[ncm] = 0.0;
    cc[ncm] = 1.0;
    ff[ncm] = 0.0;

    i = prog_right(nc,aa,bb,cc,ff,al,y1);
    if (i!=0) return i;

    for (i=0; i<ncm; i++) ff[i] = 0.0;
    ff[ncm] = 1.0;

    i = prog_right(nc,aa,bb,cc,ff,al,y2);
    if (i!=0) return i;

    ff[0] = 1.0;
    for (i=1; i<=ncm; i++) ff[i] = 0.0;

    i = prog_right(nc,aa,bb,cc,ff,al,y3);
    if (i!=0) return i;
  }
  else {
    aa[0] = 0.0;
    bb[0] = 0.0;
    cc[0] = 1.0;
    ff[0] = 0.0;

    i = prog_right(nc,aa,bb,cc,ff,al,y1);
    if (i!=0) return i;

    ff[0] = 1.0;
    for (i=1; i<=ncm; i++) ff[i] = 0.0;

    i = prog_right(nc,aa,bb,cc,ff,al,y3);
    if (i!=0) return i;
  }

  aa[0] = a0;
  bb[0] = b0;
  cc[0] = c0;
  ff[0] = f0;

  aa[ncm] = a1;
  bb[ncm] = b1;
  cc[ncm] = c1;
  ff[ncm] = f1;

  for (i=0; i<4*ncp; i++) dd[i] = 0;
  for (i=0; i<4*ncp; i++) ee[i] = 0;

  if (mp==0) {
    c1 = c1 - a1 * y2[ncm-1];
    f1 = f1 + a1 * y1[ncm-1];
    a1 = 0.0;

    dd[0] = a1;
    dd[1] = b1;
    dd[2] = c1;
    dd[3] = f1;
  }
  else if (mp<np-1) {
    c0 = c0 - b0 * y3[1];
    f0 = f0 + b0 * y1[1];
    b0 = b0 * y2[1];

    c1 = c1 - a1 * y2[ncm-1];
    f1 = f1 + a1 * y1[ncm-1];
    a1 = a1 * y3[ncm-1];

    i = mp * 8 - 4;
    dd[i]   = a0;
    dd[i+1] = b0;
    dd[i+2] = c0;
    dd[i+3] = f0;
    dd[i+4] = a1;
    dd[i+5] = b1;
    dd[i+6] = c1;
    dd[i+7] = f1;
  }
  else {
    c0 = c0 - b0 * y3[1];
    f0 = f0 + b0 * y1[1];
    b0 = 0.0;

    i = mp * 8 - 4;
    dd[i]   = a0;
    dd[i+1] = b0;
    dd[i+2] = c0;
    dd[i+3] = f0;
  }

  if (np<2) {
    for (i=0; i<ncp; i++) {
      j = 4*i;
      aa[i] = dd[j];
      bb[i] = dd[j+1];
      cc[i] = dd[j+2];
      ff[i] = dd[j+3];
    }
  }
  else {
    MPI_Allreduce(dd,ee,4*ncp,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    for (i=0; i<ncp; i++) {
      j = 4*i;
      aa[i] = ee[j];
      bb[i] = ee[j+1];
      cc[i] = ee[j+2];
      ff[i] = ee[j+3];
    }
  }

  i = prog_right(ncp,aa,bb,cc,ff,al,y4);
  if (i!=0) return i;

  if (mp==0){
    b1 = y4[0];

    for (i=0; i<nc; i++) y1[i] = y1[i] + b1 * y2[i];
  }
  else if (mp<np-1) {
    a1 = y4[2*mp-1];
    b1 = y4[2*mp];

    for (i=0; i<nc; i++) y1[i] = y1[i] + a1 * y3[i] + b1 * y2[i];
  }
  else {
    a1 = y4[2*mp-1];

    for (i=0; i<nc; i++) y1[i] = y1[i] + a1 * y3[i];
  }

  return 0;
}

// Parallel variant, MPI, modification 1:

int prog_rightpm(int np, int mp, int nc, int ip,
                 double *aa, double *bb, double *cc, double *ff,
                 double *al, double *y1, double *y2, double *y3, double *y4)
{
  int i, j, ncm=nc-1, ncp=2*np-2;
  double a0, b0, c0, f0, a1, b1, c1, f1;
  double *dd=y4+ncp, *ee=dd+4*ncp;

  a0 = aa[0];
  b0 = bb[0];
  c0 = cc[0];
  f0 = ff[0];

  a1 = aa[ncm];
  b1 = bb[ncm];
  c1 = cc[ncm];
  f1 = ff[ncm];

  if (mp==0) {
    aa[ncm] = 0.0;
    bb[ncm] = 0.0;
    cc[ncm] = 1.0;
    ff[ncm] = 0.0;

    i = prog_right(nc,aa,bb,cc,ff,al,y1);
    if (i!=0) return i;

    if (ip<2) {
      for (i=0; i<ncm; i++) ff[i] = 0.0;
      ff[ncm] = 1.0;

      i = prog_right(nc,aa,bb,cc,ff,al,y2);
      if (i!=0) return i;
    }
  }
  else if (mp<np-1) {
    aa[0] = 0.0;
    bb[0] = 0.0;
    cc[0] = 1.0;
    ff[0] = 0.0;

    aa[ncm] = 0.0;
    bb[ncm] = 0.0;
    cc[ncm] = 1.0;
    ff[ncm] = 0.0;

    i = prog_right(nc,aa,bb,cc,ff,al,y1);
    if (i!=0) return i;

    if (ip<2) {
      for (i=0; i<ncm; i++) ff[i] = 0.0;
      ff[ncm] = 1.0;

      i = prog_right(nc,aa,bb,cc,ff,al,y2);
      if (i!=0) return i;

      ff[0] = 1.0;
      for (i=1; i<=ncm; i++) ff[i] = 0.0;

      i = prog_right(nc,aa,bb,cc,ff,al,y3);
      if (i!=0) return i;
    }
  }
  else {
    aa[0] = 0.0;
    bb[0] = 0.0;
    cc[0] = 1.0;
    ff[0] = 0.0;

    i = prog_right(nc,aa,bb,cc,ff,al,y1);
    if (i!=0) return i;

    if (ip<2) {
      ff[0] = 1.0;
      for (i=1; i<=ncm; i++) ff[i] = 0.0;

      i = prog_right(nc,aa,bb,cc,ff,al,y3);
      if (i!=0) return i;
    }
  }

  aa[0] = a0;
  bb[0] = b0;
  cc[0] = c0;
  ff[0] = f0;

  aa[ncm] = a1;
  bb[ncm] = b1;
  cc[ncm] = c1;
  ff[ncm] = f1;

  for (i=0; i<4*ncp; i++) dd[i] = 0;
  for (i=0; i<4*ncp; i++) ee[i] = 0;

  if (mp==0) {
    c1 = c1 - a1 * y2[ncm-1];
    f1 = f1 + a1 * y1[ncm-1];
    a1 = 0.0;

    ee[0] = a1;
    ee[1] = b1;
    ee[2] = c1;
    ee[3] = f1;
  }
  else if (mp<np-1) {
    c0 = c0 - b0 * y3[1];
    f0 = f0 + b0 * y1[1];
    b0 = b0 * y2[1];

    c1 = c1 - a1 * y2[ncm-1];
    f1 = f1 + a1 * y1[ncm-1];
    a1 = a1 * y3[ncm-1];

    i = mp * 8 - 4;
    ee[i]   = a0;
    ee[i+1] = b0;
    ee[i+2] = c0;
    ee[i+3] = f0;
    ee[i+4] = a1;
    ee[i+5] = b1;
    ee[i+6] = c1;
    ee[i+7] = f1;
  }
  else {
    c0 = c0 - b0 * y3[1];
    f0 = f0 + b0 * y1[1];
    b0 = 0.0;

    i = mp * 8 - 4;
    ee[i]   = a0;
    ee[i+1] = b0;
    ee[i+2] = c0;
    ee[i+3] = f0;
  }

  if (np<2) {
    for(i=0; i<4*ncp; i++) dd[i] = ee[i];
  }
  else {
    MPI_Allreduce(ee,dd,4*ncp,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  }

  i = prog_rightm(ncp,dd,al,y4); if (i!=0) return i;

  if (mp==0){
    b1 = y4[0];

    for (i=0; i<nc; i++) y1[i] = y1[i] + b1 * y2[i];
  }
  else if (mp<np-1) {
    a1 = y4[2*mp-1]; b1 = y4[2*mp];

    for (i=0; i<nc; i++) y1[i] = y1[i] + a1 * y3[i] + b1 * y2[i];
  }
  else {
    a1 = y4[2*mp-1];

    for (i=0; i<nc; i++) y1[i] = y1[i] + a1 * y3[i];
  }

  return 0;
}

// Parallel variant, MPI, modification 2:

int prog_rightpn(int np, int mp, MPI_Comm cm, int nc, int ip,
                 double *aa, double *bb, double *cc, double *ff,
                 double *al, double *y1, double *y2, double *y3, double *y4)
{
  int i, j, ncm=nc-1, ncp=2*np-2;
  double a0, b0, c0, f0, a1, b1, c1, f1;
  double *dd=y4+ncp, *ee=dd+4*ncp;

  if (np<2) {
    i = prog_right(nc,aa,bb,cc,ff,al,y1);
    return i;
  }

  a0 = aa[0];
  b0 = bb[0];
  c0 = cc[0];
  f0 = ff[0];

  a1 = aa[ncm];
  b1 = bb[ncm];
  c1 = cc[ncm];
  f1 = ff[ncm];

  if (mp==0) {
    aa[ncm] = 0.0;
    bb[ncm] = 0.0;
    cc[ncm] = 1.0;
    ff[ncm] = 0.0;

    i = prog_right(nc,aa,bb,cc,ff,al,y1);
    if (i!=0) return i;

    if (ip<2) {
      for (i=0; i<ncm; i++) ff[i] = 0.0;
      ff[ncm] = 1.0;

      i = prog_right(nc,aa,bb,cc,ff,al,y2);
      if (i!=0) return i;
    }
  }
  else if (mp<np-1) {
    aa[0] = 0.0;
    bb[0] = 0.0;
    cc[0] = 1.0;
    ff[0] = 0.0;

    aa[ncm] = 0.0;
    bb[ncm] = 0.0;
    cc[ncm] = 1.0;
    ff[ncm] = 0.0;

    i = prog_right(nc,aa,bb,cc,ff,al,y1);
    if (i!=0) return i;

    if (ip<2) {
      for (i=0; i<ncm; i++) ff[i] = 0.0;
      ff[ncm] = 1.0;

      i = prog_right(nc,aa,bb,cc,ff,al,y2);
      if (i!=0) return i;

      ff[0] = 1.0;
      for (i=1; i<=ncm; i++) ff[i] = 0.0;

      i = prog_right(nc,aa,bb,cc,ff,al,y3);
      if (i!=0) return i;
    }
  }
  else {
    aa[0] = 0.0;
    bb[0] = 0.0;
    cc[0] = 1.0;
    ff[0] = 0.0;

    i = prog_right(nc,aa,bb,cc,ff,al,y1);
    if (i!=0) return i;

    if (ip<2) {
      ff[0] = 1.0;
      for (i=1; i<=ncm; i++) ff[i] = 0.0;

      i = prog_right(nc,aa,bb,cc,ff,al,y3);
      if (i!=0) return i;
    }
  }

  aa[0] = a0;
  bb[0] = b0;
  cc[0] = c0;
  ff[0] = f0;

  aa[ncm] = a1;
  bb[ncm] = b1;
  cc[ncm] = c1;
  ff[ncm] = f1;

  for (i=0; i<4*ncp; i++) dd[i] = 0;
  for (i=0; i<4*ncp; i++) ee[i] = 0;

  if (mp==0) {
    c1 = c1 - a1 * y2[ncm-1];
    f1 = f1 + a1 * y1[ncm-1];
    a1 = 0.0;

    ee[0] = a1;
    ee[1] = b1;
    ee[2] = c1;
    ee[3] = f1;
  }
  else if (mp<np-1) {
    c0 = c0 - b0 * y3[1];
    f0 = f0 + b0 * y1[1];
    b0 = b0 * y2[1];

    c1 = c1 - a1 * y2[ncm-1];
    f1 = f1 + a1 * y1[ncm-1];
    a1 = a1 * y3[ncm-1];

    i = mp * 8 - 4;
    ee[i]   = a0;
    ee[i+1] = b0;
    ee[i+2] = c0;
    ee[i+3] = f0;
    ee[i+4] = a1;
    ee[i+5] = b1;
    ee[i+6] = c1;
    ee[i+7] = f1;
  }
  else {
    c0 = c0 - b0 * y3[1];
    f0 = f0 + b0 * y1[1];
    b0 = 0.0;

    i = mp * 8 - 4;
    ee[i]   = a0;
    ee[i+1] = b0;
    ee[i+2] = c0;
    ee[i+3] = f0;
  }

  if (np<2) {
    for(i=0; i<4*ncp; i++) dd[i] = ee[i];
  }
  else {
    MPI_Allreduce(ee,dd,4*ncp,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  }

  i = prog_rightm(ncp,dd,al,y4);
  if (i!=0) return i;

  if (mp==0){
    b1 = y4[0];

    for (i=0; i<nc; i++) y1[i] = y1[i] + b1 * y2[i];
  }
  else if (mp<np-1) {
    a1 = y4[2*mp-1];
    b1 = y4[2*mp];

    for (i=0; i<nc; i++) y1[i] = y1[i] + a1 * y3[i] + b1 * y2[i];
  }
  else {
    a1 = y4[2*mp-1];

    for (i=0; i<nc; i++) y1[i] = y1[i] + a1 * y3[i];
  }

  return 0;
}