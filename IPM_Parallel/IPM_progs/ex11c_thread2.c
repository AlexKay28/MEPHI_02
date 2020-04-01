arg->np#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "mycom.h"
#include "mynet.h"
#include "myprog.h"
#include <pthread.h>
#define pi 3.14159265358979323846264338327950280

int np, mp, nl, ier, lp;
char pname[MPI_MAX_PROCESSOR_NAME];
char sname[10] = "ex11c.p00";
MPI_Status status;
union_t buf;
double tick, t1, t2, t3;

pthread_t* threads;

FILE *Fi = NULL;
FILE *Fo = NULL;

int nx;
double xa, xb, ua, ub, ak, px, px2;

double k(double x);
double k(double x) {
  //return 1.0 + ak*(x-xa)*(x-xa);
  return 1.0 + exp(-5 * ak*(x-xa));
}

double k1(double x);
double k1(double x) {
  //return ak*2.0*(x-xa);
  return -5*ak*exp(-5*ak*(x-xa));
}

double q(double x);
double q(double x) {
  //return 1.0 + ak*(xb-x)*(xb-x);
  return 1.0 - 0.5*sin(10*ak*(x-xa));
}

// THIS IS TEST FUNCTION, idk that
double u(double x);
double u(double x) {
  //return ua*cos(px*(x-xa)) + ub*sin(px*(x-xa));
  return ua*cos(px*(x-xa)) - (1/px)*ub*sin(px*(x-xa));
}

double u1(double x);
double u1(double x) {
  //return px*(-ua*sin(px*(x-xa)) + ub*cos(px*(x-xa)));
  return px*(-ua*sin(px*(x-xa)) - (1/px)*ub*cos(px*(x-xa)));
}

double u2(double x);
double u2(double x) {
  //return -px2*(ua*cos(px*(x-xa)) + ub*sin(px*(x-xa)));
  return px2*(-ua*cos(px*(x-xa)) + (1/px)*ub*sin(px*(x-xa)));
}

double f(double x);
double f(double x) {
  return -k1(x)*u1(x) - k(x)*u2(x) + q(x)*u(x);
}

typedef struct ARGS {
  int nc, ncm, ncp, np, mp;
  double *aa, *bb, *cc, *dd, *ee, *ff, *al, *y1, *y2, *y3, *y4, *xx;
  double t2;
}

void* gen_part(void *args);
void* gen_part(void *args){
  // ссылки на переменные
  someArgs_t *arg = (someArgs_t*) args;

  // случае когда всего 1 процесс..
  if (arg->arg->np<2) {
    // правая прогонка
    ier = prog_right(arg->nc,arg->aa,arg->bb,arg->cc,arg->ff,arg->al,arg->y1);
    if (ier!=0) mpierr("Bad solution 1",1);
    arg->t2 = 0.0;
  } // процессов больше чем 2
  else {
    y2 = (double*)(malloc(sizeof(double)*arg->nc));
    y3 = (double*)(malloc(sizeof(double)*arg->nc));
    y4 = (double*)(malloc(sizeof(double)*arg->ncp));
    arg->dd = (double*)(malloc(sizeof(double)*4*arg->ncp));
    arg->ee = (double*)(malloc(sizeof(double)*4*arg->ncp));

    a0 = arg->aa[0];   b0 = arg->bb[0];   c0 = arg->cc[0];   f0 = arg->ff[0];
    a1 = arg->aa[arg->ncm]; b1 = arg->bb[arg->ncm]; c1 = arg->cc[arg->ncm]; f1 = arg->ff[arg->ncm];

    // первый процесс
    if (arg->mp==0) {
      arg->aa[arg->ncm] = 0.0;
      arg->bb[arg->ncm] = 0.0;
      arg->cc[arg->ncm] = 1.0;
      arg->ff[arg->ncm] = 0.0;

      ier = prog_right(arg->nc,arg->aa,arg->bb,arg->cc,arg->ff,arg->al,arg->y1);
      if (ier!=0) mpierr("Bad solution 1",1);

      for (i=0; i<arg->ncm; i++) arg->ff[i] = 0.0;
      arg->ff[arg->ncm] = 1.0;
      ier = prog_right(arg->nc,arg->aa,arg->bb,arg->cc,arg->ff,arg->al,y2);
      if (ier!=0) mpierr("Bad solution 2",2);

      if (lp>0)
        for (i=0; i<arg->nc; i++)
          fprintf(Fo,"i=%8d x=%12le arg->y1=%12le y2=%12le\n",
            i,arg->arg->xx[i],arg->y1[i],y2[i]);
    }
    // процессы на внутренних узлах сетки
    else if (arg->mp<arg->np-1) {
      arg->aa[0] = 0.0;
      arg->bb[0] = 0.0;
      arg->cc[0] = 1.0;
      arg->ff[0] = 0.0;

      arg->aa[arg->ncm] = 0.0;
      arg->bb[arg->ncm] = 0.0;
      arg->cc[arg->ncm] = 1.0;
      arg->ff[arg->ncm] = 0.0;

      ier = prog_right(arg->nc,arg->aa,arg->bb,arg->cc,arg->ff,arg->al,arg->y1);
      if (ier!=0) mpierr("Bad solution 1",1);

      for (i=0; i<arg->ncm; i++) arg->ff[i] = 0.0;
      arg->ff[arg->ncm] = 1.0;

      ier = prog_right(arg->nc,arg->aa,arg->bb,arg->cc,arg->ff,arg->al,y2);
      if (ier!=0) mpierr("Bad solution 2",2);

      arg->ff[0] = 1.0;
      for (i=1; i<=arg->ncm; i++) arg->ff[i] = 0.0;

      ier = prog_right(arg->nc,arg->aa,arg->bb,arg->cc,arg->ff,arg->al,y3);
      if (ier!=0) mpierr("Bad solution 3",3);

      if (lp>0)
        for (i=0; i<arg->nc; i++)
          fprintf(Fo,"i=%8d x=%12le arg->y1=%12le y2=%12le y3=%12le\n",
            i,arg->xx[i],arg->y1[i],y2[i],y3[i]);
    } // процесс на последнем узле сетки
    else {
      arg->aa[0] = 0.0;
      arg->bb[0] = 0.0;
      arg->cc[0] = 1.0;
      arg->ff[0] = 0.0;

      ier = prog_right(arg->nc,arg->aa,arg->bb,arg->cc,arg->ff,arg->al,arg->y1);
      if (ier!=0) mpierr("Bad solution 1",1);

      arg->ff[0] = 1.0;
      for (i=1; i<=arg->ncm; i++) arg->ff[i] = 0.0;
      ier = prog_right(arg->nc,arg->aa,arg->bb,arg->cc,arg->ff,arg->al,y3);
      if (ier!=0) mpierr("Bad solution 3",3);

      if (lp>0)
        for (i=0; i<arg->nc; i++)
          fprintf(Fo,"i=%8d x=%12le arg->y1=%12le y3=%12le\n",
            i,arg->xx[i],arg->y1[i],y3[i]);
    }

    // инициализация-наполнение массивов arg->dd and arg->ee
    for (i=0; i<4*arg->ncp; i++) arg->dd[i] = 0;
    for (i=0; i<4*arg->ncp; i++) arg->ee[i] = 0;

    if (arg->mp==0) {
      // A0 + C0  = F0
      c1 = c1 - a1 * y2[arg->ncm-1];
      f1 = f1 + a1 * arg->y1[arg->ncm-1];
      a1 = 0.0;
      arg->dd[0] = a1;
      arg->dd[1] = b1;
      arg->dd[2] = c1;
      arg->dd[3] = f1;
    } // A + C + B = F
    else if (arg->mp<arg->np-1) {
      c0 = c0 - b0 * y3[1];
      f0 = f0 + b0 * arg->y1[1];
      b0 = b0 * y2[1];
      c1 = c1 - a1 * y2[arg->ncm-1];
      f1 = f1 + a1 * arg->y1[arg->ncm-1];
      a1 = a1 * y3[arg->ncm-1];
      i = arg->mp * 8 - 4;
      arg->dd[i]   = a0;
      arg->dd[i+1] = b0;
      arg->dd[i+2] = c0;
      arg->dd[i+3] = f0;
      arg->dd[i+4] = a1;
      arg->dd[i+5] = b1;
      arg->dd[i+6] = c1;
      arg->dd[i+7] = f1;
    }
    else {
      c0 = c0 - b0 * y3[1];
      f0 = f0 + b0 * arg->y1[1];
      b0 = 0.0;
      i = arg->mp * 8 - 4;
      arg->dd[i]   = a0;
      arg->dd[i+1] = b0;
      arg->dd[i+2] = c0;
      arg->dd[i+3] = f0;
    }

    arg->t2 = MPI_Wtime();

    // Combines values from all processes and distributes the result back to all processes
    MPI_Allreduce(arg->dd,arg->ee,4*arg->ncp,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    arg->t2 = MPI_Wtime() - arg->t2;

    if (lp>0)
      for (i=0; i<4*arg->ncp; i++)
        fprintf(Fo,"i=%8d d=%12le e=%12le\n", i, arg->dd[i], arg->ee[i]);

    for (i=0; i<arg->ncp; i++) {
      j = 4*i;
      arg->aa[i] = arg->ee[j];
      arg->bb[i] = arg->ee[j+1];
      arg->cc[i] = arg->ee[j+2];
      arg->ff[i] = arg->ee[j+3];
    }

    ier = prog_right(arg->ncp,arg->aa,arg->bb,arg->cc,arg->ff,arg->al,y4);
    if (ier!=0) mpierr("Bad solution 4",4);

    if (arg->mp==0){
      b1 = y4[0];
      for (i=0; i<arg->nc; i++)
        arg->y1[i] = arg->y1[i] + b1 * y2[i];
    }
    else if (arg->mp<arg->arg->np-1) {
      a1 = y4[2*arg->mp-1]; b1 = y4[2*arg->mp];
      for (i=0; i<arg->nc; i++)
        arg->y1[i] = arg->y1[i] + a1 * y3[i] + b1 * y2[i];
    }
    else {
      a1 = y4[2*arg->mp-1];
      for (i=0; i<arg->nc; i++)
        arg->y1[i] = arg->y1[i] + a1 * y3[i];
    }
  }
}

int main(int argc, char *argv[])
{
  int i, j, i1, i2, nc, ncm, ncp, ncx, nt;
  double hx, hx2, s0, s1, s2, a0, b0, c0, f0, a1, b1, c1, f1;
  double *xx, *aa, *bb, *cc, *dd, *ee, *ff, *al, *y1, *y2, *y3, *y4;

  someArgs_t ARGS[3];

  MyNetInit(&argc, &argv, &np, &mp, &nl, pname, &tick);


  fprintf(stderr,"Netsize: %d, process: %d, system: %s, tick=%12le\n",np,mp,pname,tick);
  sleep(1);

  sprintf(sname+7, "%02d", mp);
  ier = fopen_m(&Fo, sname, "wt");
  if (ier!=0) mpierr("Protocol file not opened",1);

  if (mp==0) {
    ier = fopen_m(&Fi,"ex11a.d","rt");
    if (ier!=0) mpierr("Data file not opened",2);
    i = fscanf(Fi,"xa=%le\n",&xa);
    i = fscanf(Fi,"xb=%le\n",&xb);
    i = fscanf(Fi,"ua=%le\n",&ua);
    i = fscanf(Fi,"ub=%le\n",&ub);
    i = fscanf(Fi,"nx=%d\n",&nx);
    i = fscanf(Fi,"lp=%d\n",&lp);
    fclose_m(&Fi);
    if (argc>1) sscanf(argv[1],"%d",&nx);
  }

  if (np>1) {
    if (mp==0) {
      buf.ddata[0] = xa; buf.ddata[1] = xb;
      buf.ddata[2] = ua; buf.ddata[3] = ub;
      buf.idata[8] = nx; buf.idata[9] = lp;
    }

    // рассылка информации от одного процесса всем остальным членам некоторой области связи
    MPI_Bcast(buf.ddata,5,MPI_DOUBLE,0,MPI_COMM_WORLD);

    if (mp>0) {
      xa = buf.ddata[0]; xb = buf.ddata[1];
      ua = buf.ddata[2]; ub = buf.ddata[3];
      nx = buf.idata[8]; lp = buf.idata[9];
    }
  }

  fprintf(Fo,"Netsize: %d, process: %d, system: %s, tick=%12le\n",np,mp,pname,tick);
  fprintf(Fo,"xa=%le xb=%le ua=%le ub=%le nx=%d lp=%d\n",xa,xb,ua,ub,nx,lp);

  t1 = MPI_Wtime();

  ak = 1.0/(xb-xa);
  // px = 0.5 pi/(xb-xa)
  px = pi/(xb-xa); // вытекает из тестовой функции
  px2 = px*px;

  hx = (xb-xa)/nx;
  hx2 = hx * hx;

  // важно получение i1, i2, nc...
  MyRange(np, mp, 0, nx, &i1, &i2, &nc);
  ncm = nc-1; // MPI??
  ncp = 2*(np-1);
  ncx = imax(nc,ncp);

  fprintf("np=%d mp=%d nx=%d ncm=%d ncp=%d ncx=%d\n",np, mp, nx, ncm, ncp, ncx);
  fprintf(Fo,"i1=%d i2=%d nc=%d\n",i1,i2,nc);

  xx = (double*)(malloc(sizeof(double)*nc));
  aa = (double*)(malloc(sizeof(double)*ncx));
  bb = (double*)(malloc(sizeof(double)*ncx));
  cc = (double*)(malloc(sizeof(double)*ncx));
  ff = (double*)(malloc(sizeof(double)*ncx));
  al = (double*)(malloc(sizeof(double)*ncx));
  y1 = (double*)(malloc(sizeof(double)*nc));

  // я думаю тут можно добавить пару тредов!
  for (i=0; i<nc; i++)
    xx[i] = xa + hx * (i1 + i); // split out the work space on grid

  // для первого узла процесса
  // если это первый процесс - следовательно первый узел удовлетворяет левому г.у
  if (mp==0) {
    aa[0] = 0.0;
    bb[0] = 0.0;
    cc[0] = 1.0;
    ff[0] = ua;
  }
  else {
    s0 = k(xx[0]);
    s1 = k(xx[0]-hx);
    s2 = k(xx[0]+hx);

    aa[0] = 0.5 * (s0 + s1);
    bb[0] = 0.5 * (s0 + s2);
    cc[0] = hx2 * q(xx[0]) + aa[0] + bb[0];
    ff[0] = hx2 * f(xx[0]);
  }

  // для внутренних узлов процесса
  // я думаю тут можно добавить пару тредов!
  for (i=1; i<ncm; i++) {
    s0 = k(xx[i]);
    s1 = k(xx[i-1]);
    s2 = k(xx[i+1]);
    aa[i] = 0.5 * (s0 + s1);
    bb[i] = 0.5 * (s0 + s2);
    cc[i] = hx2 * q(xx[i]) + aa[i] + bb[i];
    ff[i] = hx2 * f(xx[i]);
  }

  // последний узел сетки
  // если процесс последний, следовательно его последний узел удовлетворяет првым граничным условиям
  if (mp==np-1) {
    s0 = k(xx[ncm]);
    s1 = k(xx[ncm-1]);
    aa[ncm] = 0.5 * (s0 + s1);
    bb[ncm] = 0.0;
    cc[ncm] = 0.5 * hx2 * q(xx[ncm]) + aa[ncm];
    ff[ncm] = 0.5 * hx2 * f(xx[ncm]) + hx * ub * s0;
  }
  else {
    s0 = k(xx[ncm]);
    s1 = k(xx[ncm]-hx);
    s2 = k(xx[ncm]+hx);
    aa[ncm] = 0.5 * (s0 + s1);
    bb[ncm] = 0.5 * (s0 + s2);
    cc[ncm] = hx2 * q(xx[ncm]) + aa[ncm] + bb[ncm];
    ff[ncm] = hx2 * f(xx[ncm]);
  }

  if (lp>0)
    for (i=0; i<nc; i++)
      fprintf(Fo,"i=%8d a=%12le b=%12le c=%12le f=%12le\n", i, aa[i], bb[i], cc[i], ff[i]);

  // ================================ CENTRAL PART OF PROGRAM ==================================

  nt = 1;
  if (nt<2)
    //sum = integrate(f1,a,b,ni);

    int nc, ncm, ncp, np, mp;
    double *aa, *bb, *cc, *dd, *ee, *ff, *al, *y1, *y2, *y3, *y4, *xx;
    double t2;

    ARGS[0].aa = aa; ARGS[0].bb = bb; ARGS[0].ff = ff;
    ARGS[0].al = al; ARGS[0].y1 = y1; ARGS[0].y2 = y2; ARGS[0].y3 = y3;
    ARGS[0].y4 = y4; ARGS[0].dd = dd; ARGS[0].ee = ee; ARGS[0].cc = cc;

    ARGS[0].ncm = ncm; ARGS[0].mp = mp; ARGS[0].np = np;
    ARGS[0].ncp = ncp; ARGS[0].nc = nc;
    ARGS[0].t2 = t2;

    ARGS[0].xx = xx;
    ARGS[0].dd = dd;
    gen_part(ARGS[0]);
  else {
    if (!(threads = (pthread_t*) malloc(nt*sizeof(pthread_t))))
      myerr("server: not enough memory",1);

    for (i=0; i<nt; i++){
      // необходио определить нуобходимые переменные для данного участка кода
      // тут необходимо произвести разбиение участка на под-участки для тредов
      ARGS[i].var1 = var1; ARGS[i].var1 = var1;
      ARGS[i].var2 = var2; ARGS[i].var2 = var2;
      if (pthread_create(threads+i,0, gen_part, (void*) &ARGS[i])))
        myerr("server: cannot create thread",2);
    }

    for (i=0; i<nt; i++)
      if (pthread_join(threads[i],0)) myerr("server: cannot wait thread",3);

    free(threads);
  }

  // ================================ CENTRAL PART OF PROGRAM ==================================

  t1 = MPI_Wtime() - t1;

  s0 = 0.0;
  for (i=0; i<nc; i++) {
    s1 = u(xx[i]);
    s2 = dabs(s1-y1[i]);
    s0 = dmax(s0,s2);
    if (lp>0)
      fprintf(Fo,"i=%8d x=%12le y=%12le u=%12le d=%12le\n",
        i,xx[i],y1[i],s1,s2);
  }

  if (np>1) {
    s1 = s0;
    MPI_Allreduce(&s1, &s0, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  }

  if (mp==0) fprintf(stderr,"nx=%d t1=%le t2=%le dmax=%le\n",nx,t1,t2,s0);
  fprintf(Fo,"t1=%le t2=%le dmax=%le\n",t1,t2,s0);

  ier = fclose_m(&Fo);

  MPI_Finalize();
  return 0;
}
