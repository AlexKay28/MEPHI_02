int prog_right(int n,
               double* a, double* b, double* c, double* f,
               double* al, double* bt);

int prog_left(int n,
              double* a, double* b, double* c, double* f,
              double* al, double* bt);

int prog_meet(int n, int m,
              double* a, double* b, double* c, double* f,
              double* al, double* bt);

int prog_circle_right(int n,
                      double* a, double* b, double* c, double* f,
                      double* pp, double* qq, double* gm, double* al, double* bt);

int prog_integr_right(int n,
                      double* r0, double* r1, double* a, double* b, double* c, double* f,
                      double* y1, double* y2, double* y3, double* al, double* bt);

int prog_rightm(int n, double* a, double* al, double* y);

int prog_rightp(int np, int mp, int nc,
                double *aa, double *bb, double *cc, double *ff,
                double *al, double *y1, double *y2, double *y3,
                double *y4, double *dd, double *ee);

int prog_rightpm(int np, int mp, int nc, int ip,
                 double *aa, double *bb, double *cc, double *ff,
                 double *al, double *y1, double *y2, double *y3, double *y4);

int prog_rightpn(int np, int mp, MPI_Comm cm, int nc, int ip,
                 double *aa, double *bb, double *cc, double *ff,
                 double *al, double *y1, double *y2, double *y3, double *y4);
