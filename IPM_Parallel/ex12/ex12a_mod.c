//  Nonlinear boundary problem:
//
//  (k(u)u')' - q(u)u = - f(u), xa < x < xb
//
//  u(xa) = ua, u(xb) = ub
//
//  Simple iterations
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "mpi.h"
#include "mycom.h"
#include "mynet.h"
#include "myprog.h"

// Network objects:

int np, mp, nl, ier, lp;
char pname[MPI_MAX_PROCESSOR_NAME];
char sname[14] = "ex12a_mod.p00";
MPI_Status status;
union_t buf;
double tick, t1, t2, t3;

// Input/Output objects:

FILE *Fi = NULL;
FILE *Fo = NULL;

// Application objects:

int nx, it, itm;
double xa, xb, ua, ub, alf, eps;
double pa, pa2, uc, hx, hx2, rka, dm;

// Function prototypes:

double k(double u);
double k_u(double u);
double q(double u);
double f(double u);

double u(double x);
double u1(double u);
double u2(double u);

int MyJob();

// Function implementation:

double k(double u) {
  return 1.0 + u*u;
}

double k_u(double u) {
  return 2.0*u;
}

double q(double u) {
  return 1.0 / (1.0 + u*u);
}

double f(double u) {
  return q(u)*u - k_u(u)*u1(u)*u1(u) - k(u)*u2(u);
}

double u(double x) {
  return ua*exp(pa*(x-xa));
}

double u1(double u) {
  return pa*u;
}

double u2(double u) {
  return pa2*u;
}

// Main function:

int main(int argc, char *argv[])
{
  int i, myrc;
  double t1;

  MyNetInit(&argc,&argv,&np,&mp,&nl,pname,&tick);

  fprintf(stderr,"Netsize: %d, process: %d, system: %s, tick=%12le\n",np,mp,pname,tick);

  sleep(1);

  sprintf(sname+11,"%02d",mp);
  ier = fopen_m(&Fo,sname,"wt");
  if (ier!=0) mpierr("Protocol file not opened",1);

  if (mp==0) {
    ier = fopen_m(&Fi,"ex12a_mod.d","rt");
    if (ier!=0) mpierr("Data file not opened",2);
    i = fscanf(Fi,"xa=%le\n",&xa);
    i = fscanf(Fi,"xb=%le\n",&xb);
    i = fscanf(Fi,"ua=%le\n",&ua);
    i = fscanf(Fi,"ub=%le\n",&ub);
    i = fscanf(Fi,"alf=%le\n",&alf);
    i = fscanf(Fi,"eps=%le\n",&eps);
    i = fscanf(Fi,"itm=%d\n",&itm);
    i = fscanf(Fi,"nx=%d\n",&nx);
    i = fscanf(Fi,"lp=%d\n",&lp);
    fclose_m(&Fi);

    if (argc>1) sscanf(argv[1],"%d",&nx);
  }

  if (np>1) {
    if (mp==0) {
      buf.ddata[ 0] = xa;
      buf.ddata[ 1] = xb;
      buf.ddata[ 2] = ua;
      buf.ddata[ 3] = ub;
      buf.ddata[ 4] = alf;
      buf.ddata[ 5] = eps;
      buf.idata[12] = itm;
      buf.idata[13] = nx;
      buf.idata[14] = lp;
    }

    MPI_Bcast(buf.ddata,8,MPI_DOUBLE,0,MPI_COMM_WORLD);

    if (mp>0) {
      xa  = buf.ddata[ 0];
      xb  = buf.ddata[ 1];
      ua  = buf.ddata[ 2];
      ub  = buf.ddata[ 3];
      alf = buf.ddata[ 4];
      eps = buf.ddata[ 5];
      itm = buf.idata[12];
      nx  = buf.idata[13];
      lp  = buf.idata[14];
    }
  }

  ub = ua * exp(alf);

  fprintf(Fo,"Netsize: %d, process: %d, system: %s, tick=%12le\n",
    np,mp,pname,tick);

  fprintf(Fo,"xa=%le xb=%le ua=%le ub=%le alf=%le\n",xa,xb,ua,ub,alf);

  fprintf(Fo,"eps=%le itm=%d nx=%d lp=%d\n",eps,itm,nx,lp);

  t1 = MPI_Wtime();

  pa  = alf / (xb-xa);
  pa2 = pa*pa;
  uc  = (ub-ua) / (xb-xa);

  hx  = (xb-xa)/nx;
  hx2 = hx * hx;

  myrc = MyJob();

  t1 = MPI_Wtime() - t1;

  if (mp==0)
    fprintf(stderr,"np=%d nx=%d it=%d rc=%d time=%le dmax=%le\n",np,nx,it,myrc,t1,dm);

  fprintf(Fo,"np=%d nx=%d it=%d rc=%d time=%le dmax=%le\n",np,nx,it,myrc,t1,dm);

  ier = fclose_m(&Fo);

  MPI_Finalize();

  return 0;
}

// Main procedure:

int MyJob()
{
  int i, j, i1, i2, nc, ncm, ncp, ncx;

  double s0, s1, s2, y0_s, y0_l, y0_r, y0_m, y0_p;
  double *xx, *aa, *bb, *cc, *dd, *ee, *ff, *al, *y0, *y1, *y2, *y3, *y4;

  if (nx < 3*np) return -1;

  MyRange(np,mp,0,nx,&i1,&i2,&nc);

  ncm = nc-1;
  ncp = 2*(np-1);
  ncx = imax(nc,ncp);

  fprintf(Fo,"i1=%d i2=%d nc=%d\n",i1,i2,nc);

  xx = (double*)(malloc(sizeof(double)*nc));
  y0 = (double*)(malloc(sizeof(double)*nc));
  y1 = (double*)(malloc(sizeof(double)*nc));

  aa = (double*)(malloc(sizeof(double)*ncx));
  bb = (double*)(malloc(sizeof(double)*ncx));
  cc = (double*)(malloc(sizeof(double)*ncx));
  ff = (double*)(malloc(sizeof(double)*ncx));
  al = (double*)(malloc(sizeof(double)*ncx));

  if (np>1) {
    y2 = (double*)(malloc(sizeof(double)*nc));
    y3 = (double*)(malloc(sizeof(double)*nc));
    y4 = (double*)(malloc(sizeof(double)*ncp));
    dd = (double*)(malloc(sizeof(double)*4*ncp));
    ee = (double*)(malloc(sizeof(double)*4*ncp));
  }

  for (i=0; i<nc; i++) xx[i] = xa + hx * (i1 + i); // grid

  it = 0;

  for (i=0; i<nc; i++) y1[i] = ua + uc * (xx[i]-xa); // start solution

// Iterations:

  do {

    for (i=0; i<nc; i++) y0[i] = y1[i];

    y0_m = 0.0;
    y0_p = 0.0;

    if (np>1) {
      y0_l = y0[0];
      y0_r = y0[ncm];

      BndExch1D(np,mp,1,1,1,1,&y0_l,&y0_m,&y0_r,&y0_p);
    }

    if (mp==0) {
      aa[0] = 0.0;
      bb[0] = 0.0;
      cc[0] = 1.0;
      ff[0] = ua;
    }
    else {
      y0_s = y0[0];
      y0_r = y0[1];

      s0 = k(y0_s);
      s1 = k(y0_m);
      s2 = k(y0_r);

      aa[0] = 0.5 * (s0 + s1);
      bb[0] = 0.5 * (s0 + s2);
      cc[0] = hx2 * q(y0_s) + aa[0] + bb[0];
      ff[0] = hx2 * f(y0_s);
    }

    for (i=1; i<ncm; i++) {
      y0_s = y0[i];
      y0_l = y0[i-1];
      y0_r = y0[i+1];

      s0 = k(y0_s);
      s1 = k(y0_l);
      s2 = k(y0_r);

      aa[i] = 0.5 * (s0 + s1);
      bb[i] = 0.5 * (s0 + s2);
      cc[i] = hx2 * q(y0_s) + aa[i] + bb[i];
      ff[i] = hx2 * f(y0_s);
    }

    if (mp==np-1) {
      aa[ncm] = 0.0;
      bb[ncm] = 0.0;
      cc[ncm] = 1.0;
      ff[ncm] = ub;
    }
    else {
      y0_s = y0[ncm];
      y0_l = y0[ncm-1];

      s0 = k(y0_s);
      s1 = k(y0_l);
      s2 = k(y0_p);

      aa[ncm] = 0.5 * (s0 + s1);
      bb[ncm] = 0.5 * (s0 + s2);
      cc[ncm] = hx2 * q(y0_s) + aa[ncm] + bb[ncm];
      ff[ncm] = hx2 * f(y0_s);
    }

// Residual calulation:

    s0 = 0.0;

    for (i=0; i<nc; i++) {
      y0_s = y0[i];

      if (i==  0) { if (mp>   0) y0_l = y0_m; else y0_l = 0.0; } else y0_l = y0[i-1];
      if (i==ncm) { if (mp<np-1) y0_r = y0_p; else y0_r = 0.0; } else y0_r = y0[i+1];

      s1 = ff[i] + aa[i] * y0_l + bb[i] * y0_r - cc[i] * y0_s;

      s0 = dmax(s0,dabs(s1));
    }

    rka = s0;

    if (np>1) {
      MPI_Allreduce(&s0,&rka,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    }

// Residual output:

    if (lp>0) {
      if (mp==0)
        fprintf(stderr,"it=%d rka=%le\n",it,rka);

      fprintf(Fo,"\nit=%d rka=%le\n",it,rka);
    }

    if (rka<=eps) break;

// New iteration:

    it++;

    if (lp>0) {
      fprintf(Fo,"\n");
      for (i=0; i<nc; i++)
        fprintf(Fo,"i=%8d a=%12le b=%12le c=%12le f=%12le\n",
          (i1+i),aa[i],bb[i],cc[i],ff[i]);
    }

    if (np<2) {
      i = prog_right(nc,aa,bb,cc,ff,al,y1);
    }
    else {
      i = prog_rightp(np,mp,nc,aa,bb,cc,ff,al,y1,y2,y3,y4,dd,ee);
    }

    if (i != 0) return -2;

    if (lp>0) {
      fprintf(Fo,"\n");
      for (i=0; i<nc; i++)
        fprintf(Fo,"i=%8d y0=%12le y1=%12le dy=%12le\n",
          (i1+i),y0[i],y1[i],(y1[i]-y0[i]));
    }

  } while (it<=itm);

// Show solution:

  if (lp>0) {
    fprintf(Fo,"\n");
    for (i=0; i<nc; i++) {
      s1 = u(xx[i]);
      s2 = dabs(s1-y1[i]);
      fprintf(Fo,"i=%8d x=%12le y=%12le u=%12le d=%12le\n",
        (i1+i),xx[i],y1[i],s1,s2);
    }
    fprintf(Fo,"\n");
  }

// Accuracy:

  s0 = 0.0;

  for (i=0; i<nc; i++) {
    s1 = u(xx[i]);
    s2 = dabs(s1-y1[i]);
    s0 = dmax(s0,s2);
  }

  dm = s0;

  if (np>1) {
    MPI_Allreduce(&s0,&dm,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  }

  return 0;
}