//  Boundary problem:
//
//  (k(x)u')' - q(x) u = - f(x), xa < x < xb
//
//  u(xa) = ua, u(xb) = ub
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <pthread.h>
#include "mycom.h"
#include "mynet.h"
#include "myprog.h"

typedef struct tag_data_t {
  int nt, mt, rc;
} data_t;

pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;
int nt, ne=0, rc=0;

int np, mp, nl, ier, lp;
char pname[MPI_MAX_PROCESSOR_NAME];
char sname[13] = "ex11a_th.p00";
MPI_Status status;
union_t buf;
double tick, t1=0, t2=0, t3=0;

FILE *Fi = NULL;
FILE *Fo = NULL;

int nx;
double xa, xb, ua, ub, ak, px, px2, hx, hx2, dm=0;
double *dd, *ee, *y4;

double k(double x);
double k(double x) {
  return 1.0 + ak*(x-xa)*(x-xa);
}

double k1(double x);
double k1(double x) {
  return ak*2.0*(x-xa);
}

double q(double x);
double q(double x) {
  return 1.0 + ak*(xb-x)*(xb-x);
}

double u(double x);
double u(double x) {
  return ua*cos(px*(x-xa)) + ub*sin(px*(x-xa));
}

double u1(double x);
double u1(double x) {
  return px*(-ua*sin(px*(x-xa)) + ub*cos(px*(x-xa)));
}

double u2(double x);
double u2(double x) {
  return -px2*(ua*cos(px*(x-xa)) + ub*sin(px*(x-xa)));
}

double f(double x);
double f(double x) {
  return - k1(x)*u1(x) - k(x)*u2(x) + q(x)*u(x);
}

void* myjobt(void* d);
void* myjobt(void* d)
{
  int i, j, i1, i2, nc, ncm, ncp, ncx;
  double s0, s1, s2, a0, b0, c0, f0, a1, b1, c1, f1;
  double *xx, *aa, *bb, *cc, *ff, *al, *y1, *y2, *y3;

  data_t* data = (data_t *)d;
  int nt = data->nt;
  int mt = data->mt;

  int nn = np*nt;      // total thread number
  int mm = nt*mp + mt; // global thread id

  MyRange(nn,mm,0,nx,&i1,&i2,&nc);
  ncm = nc-1;
  ncp = 2*(nn-1);
  ncx = imax(nc,ncp);

  // пишем информацию о треде в файл
  if (lp>0) {
    fprintf(stderr,"mp=%d mt=%d mm=%d i1=%d i2=%d nc=%d\n",mp,mt,mm,i1,i2,nc);

    pthread_mutex_lock(&mut); // lock
    fprintf(Fo,"mp=%d mt=%d mm=%d i1=%d i2=%d nc=%d\n",mp,mt,mm,i1,i2,nc);
    pthread_mutex_unlock(&mut); // unlock
  }

  xx = (double*)(malloc(sizeof(double)*nc));
  aa = (double*)(malloc(sizeof(double)*ncx));
  bb = (double*)(malloc(sizeof(double)*ncx));
  cc = (double*)(malloc(sizeof(double)*ncx));
  ff = (double*)(malloc(sizeof(double)*ncx));
  al = (double*)(malloc(sizeof(double)*ncx));
  y1 = (double*)(malloc(sizeof(double)*nc));

  // оформляем массив сеточки
  for (i=0; i<nc; i++)
    xx[i] = xa + hx * (i1 + i);

  // первый узел сеточки
  if (mm==0) {
    aa[0] = 0.0; bb[0] = 0.0; cc[0] = 1.0; ff[0] = ua;
  }
  else {
    s0 = k(xx[0]); s1 = k(xx[0]-hx); s2 = k(xx[0]+hx);
    aa[0] = 0.5 * (s0 + s1);
    bb[0] = 0.5 * (s0 + s2);
    cc[0] = hx2 * q(xx[0]) + aa[0] + bb[0];
    ff[0] = hx2 * f(xx[0]);
  }

  for (i=1; i<ncm; i++) {
    s0 = k(xx[i]); s1 = k(xx[i-1]); s2 = k(xx[i+1]);
    aa[i] = 0.5 * (s0 + s1);
    bb[i] = 0.5 * (s0 + s2);
    cc[i] = hx2 * q(xx[i]) + aa[i] + bb[i];
    ff[i] = hx2 * f(xx[i]);
  }

  if (mm==nn-1) {
    aa[ncm] = 0.0; bb[ncm] = 0.0; cc[ncm] = 1.0; ff[ncm] = ub;
  }
  else {
    s0 = k(xx[ncm]); s1 = k(xx[ncm]-hx); s2 = k(xx[ncm]+hx);
    aa[ncm] = 0.5 * (s0 + s1);
    bb[ncm] = 0.5 * (s0 + s2);
    cc[ncm] = hx2 * q(xx[ncm]) + aa[ncm] + bb[ncm];
    ff[ncm] = hx2 * f(xx[ncm]);
  }

  if (lp>0) {
    fprintf(stderr,"[%d,%d] Coefficients -->\n",mp,mt);
    pthread_mutex_lock(&mut); // lock
    for (i=0; i<nc; i++)
      fprintf(Fo,"i=%8d a=%12le b=%12le c=%12le f=%12le\n",
        i,aa[i],bb[i],cc[i],ff[i]);
    pthread_mutex_unlock(&mut); // unlock
    fprintf(stderr,"[%d,%d] Coefficients <--\n",mp,mt);
  }

  if (nn<2) {
    fprintf(stderr,"[%d,%d] Solution -->\n",mp,mt);

    ier = prog_right(nc,aa,bb,cc,ff,al,y1);
    if (ier!=0) {
      data->rc = 1;
    }
    else {
      s0 = 0.0;
      for (i=0; i<nc; i++) {
        s1 = u(xx[i]);
        s2 = dabs(s1-y1[i]);
        s0 = dmax(s0,s2);

        if (lp>0){
          fprintf(Fo,"i=%8d x=%12le y=%12le u=%12le d=%12le\n",
          i,xx[i],y1[i],s1,s2);
	       }
	       dm = s0;
      }
    }
    fprintf(stderr,"[%d,%d] Solution <--\n",mp,mt);
    return 0;
  }
  else {
    y2 = (double*)(malloc(sizeof(double)*nc));
    y3 = (double*)(malloc(sizeof(double)*nc));

    // global arrays:
    if (mt==0) {
      y4 = (double*)calloc(  ncp,sizeof(double));
      dd = (double*)calloc(4*ncp,sizeof(double));
      ee = (double*)calloc(4*ncp,sizeof(double));
    }

    a0 = aa[0];   b0 = bb[0];   c0 = cc[0];   f0 = ff[0];
    a1 = aa[ncm]; b1 = bb[ncm]; c1 = cc[ncm]; f1 = ff[ncm];

    // первый тред первого процесса
    if (mm==0) {
      aa[ncm] = 0.0; bb[ncm] = 0.0; cc[ncm] = 1.0; ff[ncm] = 0.0;

      if (lp>0) fprintf(stderr,"[%d,%d] Calc y1\n",mp,mt);
      ier = prog_right(nc,aa,bb,cc,ff,al,y1);
      if (ier!=0) {
        data->rc = 1;
        goto Met1;
      }

      for (i=0; i<ncm; i++) ff[i] = 0.0; ff[ncm] = 1.0;

      if (lp>0) fprintf(stderr,"[%d,%d] Calc y2\n",mp,mt);
      ier = prog_right(nc,aa,bb,cc,ff,al,y2);
      if (ier!=0) {
        data->rc = 2;
        goto Met1;
      }

      if (lp>0) {
        fprintf(stderr,"[%d,%d] Basic functions:\n",mp,mt);
        pthread_mutex_lock(&mut); // lock
        for (i=0; i<nc; i++)
          fprintf(Fo,"i=%8d x=%12le y1=%12le y2=%12le\n",
            i,xx[i],y1[i],y2[i]);
        pthread_mutex_unlock(&mut); // unlock
      }
    }
    else if (mm<nn-1) {
      aa[0] = 0.0; bb[0] = 0.0; cc[0] = 1.0; ff[0] = 0.0;
      aa[ncm] = 0.0; bb[ncm] = 0.0; cc[ncm] = 1.0; ff[ncm] = 0.0;

      if (lp>0) fprintf(stderr,"[%d,%d] Calc y1\n",mp,mt);
      ier = prog_right(nc,aa,bb,cc,ff,al,y1);
      if (ier!=0) {
        data->rc = 1;
        goto Met1;
      }

      for (i=0; i<ncm; i++) ff[i] = 0.0; ff[ncm] = 1.0;

      if (lp>0) fprintf(stderr,"[%d,%d] Calc y2\n",mp,mt);
      ier = prog_right(nc,aa,bb,cc,ff,al,y2);
      if (ier!=0) {
        data->rc = 2;
        goto Met1;
      }

      ff[0] = 1.0; for (i=1; i<=ncm; i++) ff[i] = 0.0;

      if (lp>0) fprintf(stderr,"[%d,%d] Calc y3\n",mp,mt);
      ier = prog_right(nc,aa,bb,cc,ff,al,y3);
      if (ier!=0) {
        data->rc = 3;
        goto Met1;
      }

      if (lp>0) {
        fprintf(stderr,"[%d,%d] Basic functions:\n",mp,mt);
        pthread_mutex_lock(&mut); // lock
        for (i=0; i<nc; i++)
          fprintf(Fo,"i=%8d x=%12le y1=%12le y2=%12le y3=%12le\n",
            i,xx[i],y1[i],y2[i],y3[i]);
        pthread_mutex_unlock(&mut); // unlock
      }
    }
    else {
      aa[0] = 0.0; bb[0] = 0.0; cc[0] = 1.0; ff[0] = 0.0;

      if (lp>0) fprintf(stderr,"[%d,%d] Calc y1\n",mp,mt);
      ier = prog_right(nc,aa,bb,cc,ff,al,y1);
      if (ier!=0) {
        data->rc = 1;
        goto Met1;
      }

      ff[0] = 1.0; for (i=1; i<=ncm; i++) ff[i] = 0.0;

      if (lp>0) fprintf(stderr,"[%d,%d] Calc y3\n",mp,mt);
      ier = prog_right(nc,aa,bb,cc,ff,al,y3);
      if (ier!=0) {
        data->rc = 3;
        goto Met1;
      }

      if (lp>0) {
        fprintf(stderr,"[%d,%d] Basic functions:\n",mp,mt);
        pthread_mutex_lock(&mut); // lock
        for (i=0; i<nc; i++)
          fprintf(Fo,"i=%8d x=%12le y1=%12le y3=%12le\n",
            i,xx[i],y1[i],y3[i]);
        pthread_mutex_unlock(&mut); // unlock
      }
    }

  // Check of first stage:

    Met1:
    if (lp>0) fprintf(stderr,"[%d,%d] Check basic functions:\n",mp,mt);
    pthread_mutex_lock(&mut); // lock
    if (data->rc > 0) rc++;
    ne++;
    pthread_mutex_unlock(&mut); // unlock

  // Wait of all threads:

    while(ne<nt) mydelay(5);

    if (rc > 0) return 0;

    //  ???????
    // for (i=0; i<4*ncp; i++) dd[i] = 0;
    // for (i=0; i<4*ncp; i++) ee[i] = 0;

  // Coefficients of short system:

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

  // Sum for MPI-process:

    pthread_mutex_lock(&mut); // lock
    for (i=0; i<4*ncp; i++) ee[i] += dd[i];
    ne++;
    pthread_mutex_unlock(&mut); // unlock

  // Wait of all threads:

    while(ne<2*nt) mydelay(5);

  // Global data exchange:

    if ((np>1) && (mt==0)) {
      double t = MPI_Wtime();

      for (i=0; i<4*ncp; i++) {
        dd[i] = ee[i];
        ee[i] = 0.0;
      }

      MPI_Allreduce(dd,ee,4*ncp,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

      t2 += (MPI_Wtime() - t);
    }

  // Solution of short system:

    if (mt==0) {

      if (lp>0) {
        for (i=0; i<4*ncp; i++)
          fprintf(Fo,"i=%8d d=%12le e=%12le\n",i,dd[i],ee[i]);
      }

      for (i=0; i<ncp; i++) {
        j = 4*i;
        aa[i] = ee[j];   bb[i] = ee[j+1];
        cc[i] = ee[j+2]; ff[i] = ee[j+3];
      }

      ier = prog_right(ncp,aa,bb,cc,ff,al,y4);
      if (ier!=0) {
        data->rc = 4;
        rc = 1;
      }

      ne = 0;
    }

  // Wait of first thread:

    while(ne==2*nt) mydelay(5);

    if (rc > 0) return 0;

  // Final solution:

    if (mm==0){
      b1 = y4[0];
      for (i=0; i<nc; i++)
        y1[i] = y1[i] + b1 * y2[i];
    }
    else if (mm<nn-1) {
      a1 = y4[2*mm-1];
      b1 = y4[2*mm];
      for (i=0; i<nc; i++)
        y1[i] = y1[i] + a1 * y3[i] + b1 * y2[i];
    }
    else {
      a1 = y4[2*mm-1];
      for (i=0; i<nc; i++)
        y1[i] = y1[i] + a1 * y3[i];
    }
  }

  // Residual:

  s0 = 0.0;
  for (i=0; i<nc; i++) {
    s1 = u(xx[i]);
    s2 = dabs(s1-y1[i]);
    s0 = dmax(s0,s2);

    if (lp>0) {
      pthread_mutex_lock(&mut); // lock
      fprintf(Fo,"i=%8d x=%12le y=%12le u=%12le d=%12le\n",
        i,xx[i],y1[i],s1,s2);
      pthread_mutex_unlock(&mut); // unlock
    }
  }

  // Sum of residual:

  pthread_mutex_lock(&mut); // lock
  dm = dmax(dm,s0);
  ne++;
  pthread_mutex_unlock(&mut); // unlock

  // Wait of first thread:

  while(ne<nt) mydelay(5);

  // Final data exchange:

  if ((np>1) && (mt==0)) {
    double t = MPI_Wtime();

    s1 = dm;
    MPI_Allreduce(&s1,&dm,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

    t2 += (MPI_Wtime() - t);
  }

  return 0;
}

int main(int argc, char *argv[])
{
  int i;

  pthread_t* threads;
  data_t* data;

  MyNetInit(&argc,&argv,&np,&mp,&nl,pname,&tick);

  fprintf(stderr,"Netsize: %d, process: %d, system: %s, tick=%12le\n",np,mp,pname,tick);
  sleep(1);

  sprintf(sname+10,"%02d",mp);
  ier = fopen_m(&Fo,sname,"wt");
  if (ier!=0) mpierr("Protocol file not opened",1);

  if (mp==0) {
    ier = fopen_m(&Fi,"ex11a_th.d","rt");
    if (ier!=0) mpierr("Data file not opened",2);
    i = fscanf(Fi,"xa=%le\n",&xa);
    i = fscanf(Fi,"xb=%le\n",&xb);
    i = fscanf(Fi,"ua=%le\n",&ua);
    i = fscanf(Fi,"ub=%le\n",&ub);
    i = fscanf(Fi,"nx=%d\n",&nx);
    i = fscanf(Fi,"lp=%d\n",&lp);
    fclose_m(&Fi);
    if (argc>1) sscanf(argv[1],"%d",&nx);
    nt = 1;
    if (argc>2) sscanf(argv[2],"%d",&nt);
    if (nt < 1) nt = 1;
  }

  if (np>1) {
    if (mp==0) {
      buf.ddata[0] = xa; buf.ddata[1] = xb;
      buf.ddata[2] = ua; buf.ddata[3] = ub;
      buf.idata[8] = nx; buf.idata[9] = nt;
      buf.idata[10] = lp;
    }

    MPI_Bcast(buf.ddata,6,MPI_DOUBLE,0,MPI_COMM_WORLD);

    if (mp>0) {
      xa = buf.ddata[0]; xb = buf.ddata[1];
      ua = buf.ddata[2]; ub = buf.ddata[3];
      nx = buf.idata[8]; nt = buf.idata[9];
      lp = buf.idata[10];
    }
  }

  fprintf(Fo,"Netsize: %d, process: %d, threads: %d, system: %s, tick=%12le\n",np,mp,nt,pname,tick);
  fprintf(Fo,"xa=%le xb=%le ua=%le ub=%le nx=%d lp=%d\n",xa,xb,ua,ub,nx,lp);

  t1 = MPI_Wtime();

  ak = 1.0/((xb-xa)*(xb-xa));
  px = 0.5*pi/(xb-xa);
  px2 = px*px;

  hx = (xb-xa)/nx;
  hx2 = hx * hx;

  if (!(threads = (pthread_t*) malloc(nt*sizeof(pthread_t))))
    mpierr("Not enough memory",3);

  if (!(data = (data_t*) malloc(nt*sizeof(data_t))))
    mpierr("Not enough memory",4);

  for (i=0; i<nt; i++){
    (data+i)->nt=nt;
    (data+i)->mt=i;
    (data+i)->rc=0;

    if (pthread_create(threads+i,0,myjobt,(void*)(data+i)))
      mpierr("Can not create thread",5);
  }

  for (i=0; i<nt; i++){
    if (pthread_join(threads[i],0))
      mpierr("Can not wait thread",6);
    else {
      int k = (data+i)->rc;
      fprintf(stderr,"mp=%d mt=%d rc=%d rk=%d\n",mp,i,rc,k);
    }
  }

  free(data); free(threads);

  t1 = MPI_Wtime() - t1;

  t3 = t1 + t2;

  if (mp==0)
    fprintf(stderr,"nx=%d np=%d nt=%d t1=%le t2=%le t3=%le dmax=%le\n",nx,np,nt,t1,t2,t3,dm);

  fprintf(Fo,"t1=%le t2=%le t3=%le dmax=%le\n",t1,t2,t3,dm);

  ier = fclose_m(&Fo);

  MPI_Finalize();

  return 0;
}
