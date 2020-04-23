//  Nonlinear boundary problem:
//
//  (k(u)u')' - q(u)u = - f(u), xa < x < xb
//
//  u(xa) = ua, u(xb) = ub
//
//  Simple iterations, MPI + PThreads
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <pthread.h>
#include <mpi.h>
#include "mycom.h"
#include "mynet.h"
#include "myprog.h"

// Thread objects:
typedef struct tag_data_t {
  int nt, mt, rc;
} data_t;

int nth=1, nth_max=1;
pthread_mutex_t mut;
pthread_barrier_t bstart;

// Network objects:

int np, mp, nl, ier, lp;
char pname[MPI_MAX_PROCESSOR_NAME];
char sname[14] = "ex12a_omp.p00";
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
double *dd, *ee, *y4, *y5; 

// Function prototypes:

double k(double u);
double k_u(double u);
double q(double u);
double f(double u);

double u(double x);
double u1(double u);
double u2(double u);

void* MyJob(void* d);

int prog_right_mpi_pth(int np, int mp, int nt, int mt, int nc,
                       double *aa, double *bb, double *cc, double *ff,
                       double *al, double *y1, double *y2, double *y3,
                       double *y4, double *dd, double *ee);

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
  int i, myrc=0;
  double t1;

  pthread_t* threads;
  data_t* data;

  MyNetInit(&argc,&argv,&np,&mp,&nl,pname,&tick);

  nth_max = sysconf(_SC_NPROCESSORS_ONLN);

  fprintf(stderr,"Netsize: %d, process: %d, system: %s, max_thr: %d, tick=%12le\n",
    np, mp, pname, nth_max, tick);

  sleep(1);

  sprintf(sname+11,"%02d",mp);
  ier = fopen_m(&Fo,sname,"wt");
  if (ier!=0) mpierr("Protocol file not opened",1);

  if (mp==0) {
    ier = fopen_m(&Fi,"ex12a_omp.d","rt");
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
    if (argc>2) sscanf(argv[2],"%d",&nth);
    if (argc>3) sscanf(argv[3],"%d",&lp);

    if (nth<1) nth = 1;
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
      buf.idata[15] = nth;
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
      nth = buf.idata[15];
    }
  }

  ub = ua * exp(alf);

  fprintf(Fo,"Netsize: %d, process: %d, system: %s, max_thr: %d, thr: %d, tick=%12le\n",
    np, mp, pname, nth_max, nth, tick);

  fprintf(Fo,"xa=%le xb=%le ua=%le ub=%le alf=%le\n",xa,xb,ua,ub,alf);

  fprintf(Fo,"eps=%le itm=%d nx=%d lp=%d\n",eps,itm,nx,lp);

// Start main calculations:

  t1 = MPI_Wtime();

  pa  = alf / (xb-xa);
  pa2 = pa*pa;
  uc  = (ub-ua) / (xb-xa);

  hx  = (xb-xa)/nx;
  hx2 = hx * hx;

// Parallelizm init:

  if (!(threads = (pthread_t*) malloc(nth*sizeof(pthread_t))))
    mpierr("Not enough memory for threads",3);

  if (!(data = (data_t*) malloc(nth*sizeof(data_t))))
    mpierr("Not enough memory for data",4);

  if (pthread_mutex_init(&mut, NULL) != 0)
    mpierr("Mutex not inited",5);

  if (pthread_barrier_init(&bstart, NULL, nth) != 0)
    mpierr("Barrier not inited",6);

// Parallel:

  for (i=0; i<nth; i++){
    (data+i)->nt=nth;
    (data+i)->mt=i;
    (data+i)->rc=0;

    if (pthread_create(threads+i,0,MyJob,(void*)(data+i)))
      mpierr("Can not create thread",7);
  }

  for (i=0; i<nth; i++){
    if (pthread_join(threads[i],0))
      mpierr("Can not wait thread",8);
    else {
      int mt = (data+i)->mt;
      int rc = (data+i)->rc;
      if (rc!=0) myrc = rc;
      if (lp>0) fprintf(stderr,"[%d,%d] rc=%d\n", mp, mt, rc);
    }
  }

// Parallelizm destroy:

  pthread_mutex_destroy(&mut);
  free(data); free(threads);

  t1 = MPI_Wtime() - t1;

// Finish main calculations:

  {
    int nn = np*nth;

    if (mp==0) fprintf(stderr,"nn=%d np=%d nt=%d nx=%d it=%d rc=%d time=%le dmax=%le\n",nn,np,nth,nx,it,myrc,t1,dm);

    fprintf(Fo,"\nnn=%d np=%d nt=%d nx=%d it=%d rc=%d time=%le dmax=%le\n",nn,np,nth,nx,it,myrc,t1,dm);
  }

  ier = fclose_m(&Fo);

  MPI_Finalize();

  return 0;
}

// Main procedure:

void* MyJob(void* d)
{
  int i, j, i1, i2, nn, mm, nc, ncm, ncp, ncx;
  double s0, s1, s2, s3, s4, y0_s, y0_l, y0_r, y0_m, y0_p;
  double *xx, *aa, *bb, *cc, *ff, *al, *y0, *y1, *y2, *y3;

  data_t* data = (data_t *)d; 
  int nt = data->nt;
  int mt = data->mt;

// Initialization:

  nn = np*nt;
  mm = nt*mp + mt;

  if (nx < 3*nn) {
    data->rc = -1;
    return 0;
  }

  MyRange(nn,mm,0,nx,&i1,&i2,&nc);

  ncm = nc-1;
  ncp = 2*(nn-1);
  ncx = imax(nc,ncp);

// First output:

  if (lp>0) fprintf(stderr,"[%d,%d] np=%d nt=%d nn=%d mm=%d i1=%d i2=%d nc=%d\n",mp,mt,np,nt,nn,mm,i1,i2,nc);

  pthread_mutex_lock(&mut); // lock
  fprintf(Fo,"[%d,%d] np=%d nt=%d nn=%d mm=%d i1=%d i2=%d nc=%d\n",mp,mt,np,nt,nn,mm,i1,i2,nc);
  pthread_mutex_unlock(&mut); // unlock

// Memory allocation:

  xx = (double*)(malloc(sizeof(double)*nc));
  y0 = (double*)(malloc(sizeof(double)*nc));
  y1 = (double*)(malloc(sizeof(double)*nc));

  aa = (double*)(malloc(sizeof(double)*ncx));
  bb = (double*)(malloc(sizeof(double)*ncx));
  cc = (double*)(malloc(sizeof(double)*ncx));
  ff = (double*)(malloc(sizeof(double)*ncx));
  al = (double*)(malloc(sizeof(double)*ncx));

  if (nn>1) {
    y2 = (double*)(malloc(sizeof(double)*nc));
    y3 = (double*)(malloc(sizeof(double)*nc));

    if (mt==0) {
      y4 = (double*)(malloc(sizeof(double)*ncp));
      dd = (double*)(malloc(sizeof(double)*4*ncp));
      ee = (double*)(malloc(sizeof(double)*4*ncp));
      y5 = (double*)(malloc(sizeof(double)*(2*nt+2)));
    }
  }

// Init grid & start profile:

  for (i=0; i<nc; i++) xx[i] = xa + hx * (i1 + i); // grid

  for (i=0; i<nc; i++) y1[i] = ua + uc * (xx[i]-xa); // start solution

// Iterations:

  if (mt==0) it = 0;

  do {

// Reload iteration:

    for (i=0; i<nc; i++) y0[i] = y1[i];

// Neighbors data:

    y0_m = 0.0;
    y0_p = 0.0;

    if (nn>1) {
      pthread_barrier_wait(&bstart);
      if (mt==0) for (i=0; i<2*nt+2; i++) y5[i] = 0.0;
      pthread_barrier_wait(&bstart);

      y5[2*mt+1] = y0[0];
      y5[2*mt+2] = y0[ncm];

      pthread_barrier_wait(&bstart);

      if (mt==0) {
        s2 = 0.0;
        s4 = 0.0;

        if (np>1) {
          s1 = y5[1];
          s3 = y5[2*nt];
          BndExch1D(np,mp,1,1,1,1,&s1,&s2,&s3,&s4);
        }

        y5[0]      = s2;
        y5[2*nt+1] = s4;
      }

      pthread_barrier_wait(&bstart);

      y0_m = y5[2*mt];
      y0_p = y5[2*mt+3];
    }

// Main coefficients:

    if (mm==0) {
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

    if (mm==nn-1) {
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

    if (mt == 0) rka = 0.0;

    s0 = 0.0;

    for (i=0; i<nc; i++) {
      y0_s = y0[i];

      if (i==  0) y0_l = y0_m; else y0_l = y0[i-1];
      if (i==ncm) y0_r = y0_p; else y0_r = y0[i+1];

      s1 = ff[i] + aa[i] * y0_l + bb[i] * y0_r - cc[i] * y0_s;

      s0 = dmax(s0,dabs(s1));
    }

// Sborka:

    if (nn<2) {
      rka = s0;
    }
    else {
      pthread_barrier_wait(&bstart);
      pthread_mutex_lock(&mut); // lock
      if (rka<s0) rka = s0;
      pthread_mutex_unlock(&mut); // unlock
      pthread_barrier_wait(&bstart);

      if (mt==0) {
        if (np>1) {
          s0 = rka;
          MPI_Allreduce(&s0,&rka,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        }
      }
      pthread_barrier_wait(&bstart);
    }

// Residual output:

    if (lp>0) {
      fprintf(stderr,"[%d,%d] it=%d rka=%le\n",mp,mt,it,rka);
      if (mt==0) fprintf(Fo,"\nit=%d rka=%le\n",it,rka);
    }

    if (nn>1) pthread_barrier_wait(&bstart);

    if (rka<=eps) break;

// New iteration:

    if (mt==0) it++;

// Output coefficients:

    if (lp>1) {
      pthread_mutex_lock(&mut); // lock
      fprintf(Fo,"\n");
      for (i=0; i<nc; i++)
        fprintf(Fo,"i=%8d a=%12le b=%12le c=%12le f=%12le\n",
          (i1+i),aa[i],bb[i],cc[i],ff[i]);
      pthread_mutex_unlock(&mut); // unlock
    }

// Progonka:

    if (nn<2) {
      i = prog_right(nc,aa,bb,cc,ff,al,y1);
    }
    else {
      pthread_barrier_wait(&bstart);

      if (nt<2) {
        i = prog_rightp(np,mp,nc,aa,bb,cc,ff,al,y1,y2,y3,y4,dd,ee);
      }
      else {
        i = prog_right_mpi_pth(np,mp,nt,mt,nc,aa,bb,cc,ff,al,y1,y2,y3,y4,dd,ee);
      }

      pthread_barrier_wait(&bstart);
    }

    if (i != 0) {
      data->rc = -2;
      return 0;
    }

// Output new solution:

    if (lp>1) {
      pthread_mutex_lock(&mut); // lock
      fprintf(Fo,"\n");
      for (i=0; i<nc; i++)
        fprintf(Fo,"i=%8d y0=%12le y1=%12le dy=%12le\n",
          (i1+i),y0[i],y1[i],(y1[i]-y0[i]));
      pthread_mutex_unlock(&mut); // unlock
    }

// Exit from iteration loop:

  } while (it<=itm);

// Output solution:

  if (lp>1) {
    pthread_mutex_lock(&mut); // lock
    fprintf(Fo,"\n");
    for (i=0; i<nc; i++) {
      s1 = u(xx[i]);
      s2 = dabs(s1-y1[i]);
      fprintf(Fo,"i=%8d x=%12le y=%12le u=%12le d=%12le\n",
        (i1+i),xx[i],y1[i],s1,s2);
    }
    fprintf(Fo,"\n");
    pthread_mutex_unlock(&mut); // unlock
  }

// Accuracy:

  if (mt == 0) dm = 0.0;

  s0 = 0.0;

  for (i=0; i<nc; i++) {
    s1 = u(xx[i]);
    s2 = dabs(s1-y1[i]);
    s0 = dmax(s0,s2);
  }

  if (nn<2) {
    dm = s0;
  }
  else {
    pthread_barrier_wait(&bstart);
    pthread_mutex_lock(&mut); // lock
    if (dm<s0) dm = s0;
    pthread_mutex_unlock(&mut); // unlock
    pthread_barrier_wait(&bstart);

    if (mt==0) {
      if (np>1) {
        s0 = dm;
        MPI_Allreduce(&s0,&dm,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
      }
    }
    pthread_barrier_wait(&bstart);
  }

// Finish:

  data->rc = 0;
  return 0;
}

// Parallel variant, MPI + PThreads:

int prog_right_mpi_pth(int np, int mp, int nt, int mt, int nc,
                       double *aa, double *bb, double *cc, double *ff,
                       double *al, double *y1, double *y2, double *y3,
                       double *y4, double *dd, double *ee)
{
  int i, j, nn, mm, ncm, ncp;
  double a0, b0, c0, f0, a1, b1, c1, f1;

// Initialization:

  nn = np*nt;
  mm = nt*mp + mt;

  ncm = nc-1;
  ncp = 2*nn-2;

  if (mt==0) {
    for (i=0; i<4*ncp; i++) dd[i] = 0;
    for (i=0; i<4*ncp; i++) ee[i] = 0;
  }

// Save boundary coefficients:

  a0 = aa[0];
  b0 = bb[0];
  c0 = cc[0];
  f0 = ff[0];

  a1 = aa[ncm];
  b1 = bb[ncm];
  c1 = cc[ncm];
  f1 = ff[ncm];

// Calculation of base functions:

  if (mm==0) {
    aa[ncm] = 0.0;
    bb[ncm] = 0.0;
    cc[ncm] = 1.0;
    ff[ncm] = 0.0;

    i = prog_right(nc,aa,bb,cc,ff,al,y1);
    if (i!=0) {
      fprintf(stderr,"[%d,%d] y1: rc=%d\n",mp,mt,i);
      return i;
    }

    for (i=0; i<ncm; i++) ff[i] = 0.0;
    ff[ncm] = 1.0;

    i = prog_right(nc,aa,bb,cc,ff,al,y2);
    if (i!=0) {
      fprintf(stderr,"[%d,%d] y2: rc=%d\n",mp,mt,i);
      return i;
    }
  }
  else if (mm<nn-1) {
    aa[0] = 0.0;
    bb[0] = 0.0;
    cc[0] = 1.0;
    ff[0] = 0.0;

    aa[ncm] = 0.0;
    bb[ncm] = 0.0;
    cc[ncm] = 1.0;
    ff[ncm] = 0.0;

    i = prog_right(nc,aa,bb,cc,ff,al,y1);
    if (i!=0) {
      fprintf(stderr,"[%d,%d] y1: rc=%d\n",mp,mt,i);
      return i;
    }

    for (i=0; i<ncm; i++) ff[i] = 0.0;
    ff[ncm] = 1.0;

    i = prog_right(nc,aa,bb,cc,ff,al,y2);
    if (i!=0) {
      fprintf(stderr,"[%d,%d] y2: rc=%d\n",mp,mt,i);
      return i;
    }

    ff[0] = 1.0;
    for (i=1; i<=ncm; i++) ff[i] = 0.0;

    i = prog_right(nc,aa,bb,cc,ff,al,y3);
    if (i!=0) {
      fprintf(stderr,"[%d,%d] y3: rc=%d\n",mp,mt,i);
      return i;
    }
  }
  else {
    aa[0] = 0.0;
    bb[0] = 0.0;
    cc[0] = 1.0;
    ff[0] = 0.0;

    i = prog_right(nc,aa,bb,cc,ff,al,y1);
    if (i!=0) {
      fprintf(stderr,"[%d,%d] y1: rc=%d\n",mp,mt,i);
      return i;
    }

    ff[0] = 1.0;
    for (i=1; i<=ncm; i++) ff[i] = 0.0;

    i = prog_right(nc,aa,bb,cc,ff,al,y3);
    if (i!=0) {
      fprintf(stderr,"[%d,%d] y3: rc=%d\n",mp,mt,i);
      return i;
    }
  }

// Load boundary coefficients:

  aa[0] = a0;
  bb[0] = b0;
  cc[0] = c0;
  ff[0] = f0;

  aa[ncm] = a1;
  bb[ncm] = b1;
  cc[ncm] = c1;
  ff[ncm] = f1;

// Formation of short system:

  pthread_barrier_wait(&bstart);

  if (mm==0) {
    c1 = c1 - a1 * y2[ncm-1];
    f1 = f1 + a1 * y1[ncm-1];
    a1 = 0.0;

    dd[0] = a1;
    dd[1] = b1;
    dd[2] = c1;
    dd[3] = f1;
  }
  else if (mm<nn-1) {
    c0 = c0 - b0 * y3[1];
    f0 = f0 + b0 * y1[1];
    b0 = b0 * y2[1];

    c1 = c1 - a1 * y2[ncm-1];
    f1 = f1 + a1 * y1[ncm-1];
    a1 = a1 * y3[ncm-1];

    i = mm * 8 - 4;
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

    i = mm * 8 - 4;
    dd[i]   = a0;
    dd[i+1] = b0;
    dd[i+2] = c0;
    dd[i+3] = f0;
  }

  pthread_barrier_wait(&bstart);

  if (mt == 0) {
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
    if (i!=0) {
      fprintf(stderr,"[%d,%d] y4: rc=%d\n",mp,mt,i);
      return i;
    }
  }

  pthread_barrier_wait(&bstart);

  if (mm==0){
    b1 = y4[0];

    for (i=0; i<nc; i++) y1[i] = y1[i] + b1 * y2[i];
  }
  else if (mm<nn-1) {
    a1 = y4[2*mm-1]; b1 = y4[2*mm];

    for (i=0; i<nc; i++) y1[i] = y1[i] + a1 * y3[i] + b1 * y2[i];
  }
  else {
    a1 = y4[2*mm-1];

    for (i=0; i<nc; i++) y1[i] = y1[i] + a1 * y3[i];
  }

  pthread_barrier_wait(&bstart);

  return 0;
}
