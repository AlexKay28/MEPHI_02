#include <stdio.h>
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

pthread_t* thread1;
pthread_t* thread2;

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

typedef struct progonka_args {
    int nc;
    double *a, *b, *c, *f, *al, *y;
} someArgs_t;

void* progonka(void *args);
void* progonka(void *args){
  someArgs_t *arg = (someArgs_t*) args;
  ier = prog_right(arg->nc,arg->a,arg->b,arg->c,arg->f,arg->al,arg->y);
  //return ier;
  return 0;
}

int main(int argc, char *argv[])
{
  int i, j, i1, i2, nc, ncm, ncp, ncx;
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
  ncm = nc-1;
  ncp = 2*(np-1);
  ncx = imax(nc,ncp);

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

  // случае когда всего 1 процесс..
  if (np<2) {
    // правая прогонка
    ier = prog_right(nc,aa,bb,cc,ff,al,y1);
    if (ier!=0) mpierr("Bad solution 1",1);
    t2 = 0.0;
  } // процессов больше чем 2
  else {
    y2 = (double*)(malloc(sizeof(double)*nc));
    y3 = (double*)(malloc(sizeof(double)*nc));
    y4 = (double*)(malloc(sizeof(double)*ncp));
    dd = (double*)(malloc(sizeof(double)*4*ncp));
    ee = (double*)(malloc(sizeof(double)*4*ncp));

    a0 = aa[0];   b0 = bb[0];   c0 = cc[0];   f0 = ff[0];
    a1 = aa[ncm]; b1 = bb[ncm]; c1 = cc[ncm]; f1 = ff[ncm];

    // идея распараллеливания по тредам сконцентрирована на распараллеливании выполняемых функций программы
    // я решил что будет выгодно распараллелить методы прогонки используемые в программе, т. к
    // периодически они выставляются в стек, поэтому выполенние каждой прогонки в одельном треде может
    // сократить врея выполнения отдельных частей программы втрое

    // первый процесс
    if (mp==0) {
      aa[ncm] = 0.0;
      bb[ncm] = 0.0;
      cc[ncm] = 1.0;
      ff[ncm] = 0.0;
      // ==================================================== Pthread ->

      //thread 1
      ARGS[1].nc = nc; ARGS[1].a = aa; ARGS[1].b = bb;
      ARGS[1].c = cc; ARGS[1].f = ff; ARGS[1].al = al;
      ARGS[1].y = y1;
      if (pthread_create(thread1, 0, progonka, (void*) &ARGS[1])) {
        myerr("server: cannot create thread",2);
      }
      if (ier!=0) mpierr("Bad solution 1",1);

      //thread 2
      for (i=0; i<ncm; i++) ff[i] = 0.0;
      ff[ncm] = 1.0;
      ier = prog_right(nc,aa,bb,cc,ff,al,y2);
      if (ier!=0) mpierr("Bad solution 2",2);

      // pthread_join(thread,0);
      pthread_join(*thread1,0);
      free(thread1);
      // ==================================================== Pthread <-
      if (lp>0)
        for (i=0; i<nc; i++)
          fprintf(Fo,"i=%8d x=%12le y1=%12le y2=%12le\n",
            i,xx[i],y1[i],y2[i]);
    }
    // процессы на внутренних узлах сетки
    else if (mp<np-1) {
      aa[0] = 0.0;
      bb[0] = 0.0;
      cc[0] = 1.0;
      ff[0] = 0.0;

      aa[ncm] = 0.0;
      bb[ncm] = 0.0;
      cc[ncm] = 1.0;
      ff[ncm] = 0.0;

      // ==================================================== Pthread ->
      //ier = prog_right(nc,aa,bb,cc,ff,al,y1);
      // thread 1
      ARGS[1].nc = nc; ARGS[1].a = aa; ARGS[1].b = bb;
      ARGS[1].c = cc; ARGS[1].f = ff; ARGS[1].al = al;
      ARGS[1].y = y1;
      if (pthread_create(thread1, 0, progonka, (void*) &ARGS[1])) {
        myerr("server: cannot create thread",2);
      }
      if (ier!=0) mpierr("Bad solution 1",1);



      for (i=0; i<ncm; i++) ff[i] = 0.0;
      ff[ncm] = 1.0;

      //ier = prog_right(nc,aa,bb,cc,ff,al,y2);
      // thread 2
      ARGS[2].nc = nc; ARGS[2].a = aa; ARGS[2].b = bb;
      ARGS[2].c = cc; ARGS[2].f = ff; ARGS[2].al = al;
      ARGS[2].y = y2;
      if (pthread_create(thread2, 0, progonka, (void*) &ARGS[2])) {
        myerr("server: cannot create thread",2);
      }
      if (ier!=0) mpierr("Bad solution 2",2);

      ff[0] = 1.0;
      for (i=1; i<=ncm; i++) ff[i] = 0.0;

      ier = prog_right(nc,aa,bb,cc,ff,al,y3); // thread 3
      if (ier!=0) mpierr("Bad solution 3",3);

      // join thread 1
      // join thread 2
      // join thread 3
      pthread_join(*thread1,0);
      pthread_join(*thread2,0);
      free(thread1);
      free(thread2);
      // ==================================================== Pthread <-
      if (lp>0)
        for (i=0; i<nc; i++)
          fprintf(Fo,"i=%8d x=%12le y1=%12le y2=%12le y3=%12le\n",
            i,xx[i],y1[i],y2[i],y3[i]);
    } // процесс на последнем узле сетки
    else {
      aa[0] = 0.0;
      bb[0] = 0.0;
      cc[0] = 1.0;
      ff[0] = 0.0;
      // ==================================================== Pthread ->
      //ier = prog_right(nc,aa,bb,cc,ff,al,y1);
      //thread 1
      ARGS[1].nc = nc; ARGS[1].a = aa; ARGS[1].b = bb;
      ARGS[1].c = cc; ARGS[1].f = ff; ARGS[1].al = al;
      ARGS[1].y = y1;
      if (pthread_create(thread1, 0, progonka, (void*) &ARGS[1])) {
        myerr("server: cannot create thread",2);
      }
      if (ier!=0) mpierr("Bad solution 1",1);;

      ff[0] = 1.0;
      for (i=1; i<=ncm; i++) ff[i] = 0.0;

      ier = prog_right(nc,aa,bb,cc,ff,al,y3); //thread 2
      if (ier!=0) mpierr("Bad solution 3",3);
      // join thread 1
      // join thread 2
      pthread_join(*thread1,0);
      free(thread1);
      // ==================================================== Pthread <-
      if (lp>0)
        for (i=0; i<nc; i++)
          fprintf(Fo,"i=%8d x=%12le y1=%12le y3=%12le\n",
            i,xx[i],y1[i],y3[i]);
    }
    // инициализация-наполнение массивов dd and ee
    for (i=0; i<4*ncp; i++) dd[i] = 0;
    for (i=0; i<4*ncp; i++) ee[i] = 0;

    if (mp==0) {
      // A0 + C0  = F0
      c1 = c1 - a1 * y2[ncm-1];
      f1 = f1 + a1 * y1[ncm-1];
      a1 = 0.0;
      dd[0] = a1;
      dd[1] = b1;
      dd[2] = c1;
      dd[3] = f1;
    } // A + C + B = F
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

    t2 = MPI_Wtime();

    // Combines values from all processes and distributes the result back to all processes
    MPI_Allreduce(dd,ee,4*ncp,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    t2 = MPI_Wtime() - t2;

    if (lp>0)
      for (i=0; i<4*ncp; i++)
        fprintf(Fo,"i=%8d d=%12le e=%12le\n", i, dd[i], ee[i]);

    for (i=0; i<ncp; i++) {
      j = 4*i;
      aa[i] = ee[j];
      bb[i] = ee[j+1];
      cc[i] = ee[j+2];
      ff[i] = ee[j+3];
    }

    ier = prog_right(ncp,aa,bb,cc,ff,al,y4);
    if (ier!=0) mpierr("Bad solution 4",4);

    if (mp==0){
      b1 = y4[0];
      for (i=0; i<nc; i++)
        y1[i] = y1[i] + b1 * y2[i];
    }
    else if (mp<np-1) {
      a1 = y4[2*mp-1]; b1 = y4[2*mp];
      for (i=0; i<nc; i++)
        y1[i] = y1[i] + a1 * y3[i] + b1 * y2[i];
    }
    else {
      a1 = y4[2*mp-1];
      for (i=0; i<nc; i++)
        y1[i] = y1[i] + a1 * y3[i];
    }
  }

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
