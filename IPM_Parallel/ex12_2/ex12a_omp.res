IMM10:
---
>mpicc -o ex12a_omp.px -O2 -fopenmp ex12a_omp.c mycom.c mynet.c myprog.c -lm
---
>ex12a_omp.px 10..1000000 1
nn=1 np=1 nt=1 nx=10      it=30 rc=0 time=4.911423e-05 dmax=2.830711e-03
nn=1 np=1 nt=1 nx=100     it=24 rc=0 time=7.729232e-05 dmax=2.869907e-05
nn=1 np=1 nt=1 nx=1000    it=20 rc=0 time=1.282059e-03 dmax=2.871838e-07
nn=1 np=1 nt=1 nx=10000   it=16 rc=0 time=1.059446e-02 dmax=1.578371e-08
nn=1 np=1 nt=1 nx=100000  it=12 rc=0 time=4.027491e-02 dmax=9.606953e-07
nn=1 np=1 nt=1 nx=1000000 it=7  rc=0 time=2.308235e-01 dmax=2.028484e-04
---
>export OMP_NUM_THREADS=32
>ex12a_omp.px 100000 1..32
nn=1  np=1  nt=1  nx=100000 it=12 rc=0 time=4.943637e-02 dmax=9.606953e-07
nn=2  np=1  nt=2  nx=100000 it=12 rc=0 time=3.936188e-02 dmax=9.604200e-07
nn=4  np=1  nt=4  nx=100000 it=12 rc=0 time=2.414690e-02 dmax=9.596533e-07
nn=8  np=1  nt=8  nx=100000 it=12 rc=0 time=1.553845e-02 dmax=9.591225e-07
nn=16 np=1  nt=16 nx=100000 it=12 rc=0 time=1.154919e-02 dmax=9.591722e-07
nn=32 np=1  nt=32 nx=100000 it=12 rc=0 time=1.349404e-02 dmax=9.574738e-07
---
>mpirun -np 1..32 ex12a_omp.px 100000 1
nn=1  np=1  nt=1  nx=100000 it=12 rc=0 time=3.840055e-02 dmax=9.606953e-07
nn=2  np=2  nt=1  nx=100000 it=12 rc=0 time=4.043683e-02 dmax=9.604200e-07
nn=4  np=4  nt=1  nx=100000 it=12 rc=0 time=2.068361e-02 dmax=9.596533e-07
nn=8  np=8  nt=1  nx=100000 it=12 rc=0 time=2.006691e-02 dmax=9.591225e-07
nn=16 np=16 nt=1  nx=100000 it=12 rc=0 time=1.172411e-02 dmax=9.591722e-07
nn=32 np=32 nt=1  nx=100000 it=12 rc=0 time=9.085771e-03 dmax=9.574738e-07
---
>mpirun -np 1..32 ex12a_omp.px 100000 1..32
nn=1  np=1  nt=1  nx=100000 it=12 rc=0 time=5.069185e-02 dmax=9.606953e-07
-
nn=2  np=1  nt=2  nx=100000 it=12 rc=0 time=3.563725e-02 dmax=9.604200e-07
nn=2  np=2  nt=1  nx=100000 it=12 rc=0 time=3.459290e-02 dmax=9.604200e-07 ***
-
nn=4  np=1  nt=4  nx=100000 it=12 rc=0 time=4.056900e-02 dmax=9.596533e-07
nn=4  np=2  nt=2  nx=100000 it=12 rc=0 time=2.801352e-02 dmax=9.596533e-07 ***
nn=4  np=4  nt=1  nx=100000 it=12 rc=0 time=2.804458e-02 dmax=9.596533e-07
-
nn=8  np=1  nt=8  nx=100000 it=12 rc=0 time=4.229580e-02 dmax=9.591225e-07
nn=8  np=2  nt=4  nx=100000 it=12 rc=0 time=2.296102e-02 dmax=9.591225e-07
nn=8  np=4  nt=2  nx=100000 it=12 rc=0 time=1.664455e-02 dmax=9.591225e-07 ***
nn=8  np=8  nt=1  nx=100000 it=12 rc=0 time=1.865252e-02 dmax=9.591225e-07
-
nn=16 np=1  nt=16 nx=100000 it=12 rc=0 time=5.138818e-02 dmax=9.591722e-07
nn=16 np=2  nt=8  nx=100000 it=12 rc=0 time=2.494914e-02 dmax=9.591722e-07
nn=16 np=4  nt=4  nx=100000 it=12 rc=0 time=1.302756e-02 dmax=9.591722e-07
nn=16 np=8  nt=2  nx=100000 it=12 rc=0 time=1.220790e-02 dmax=9.591722e-07 ***
nn=16 np=16 nt=1  nx=100000 it=12 rc=0 time=1.225849e-02 dmax=9.591722e-07
-
nn=32 np=1  nt=32 nx=100000 it=12 rc=0 time=6.202360e-02 dmax=9.574738e-07
nn=32 np=2  nt=16 nx=100000 it=12 rc=0 time=3.072023e-02 dmax=9.574738e-07
nn=32 np=4  nt=8  nx=100000 it=12 rc=0 time=7.251263e-03 dmax=9.574738e-07 ***
nn=32 np=8  nt=4  nx=100000 it=12 rc=0 time=7.687859e-03 dmax=9.574738e-07
nn=32 np=16 nt=2  nx=100000 it=12 rc=0 time=8.718215e-03 dmax=9.574738e-07
nn=32 np=32 nt=1  nx=100000 it=12 rc=0 time=1.046699e-02 dmax=9.574738e-07
-