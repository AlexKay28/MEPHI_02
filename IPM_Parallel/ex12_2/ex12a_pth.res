IMM10:
---
>mpicc -o ex12a_pth.px -O2 -pthread ex12a_pth.c mycom.c mynet.c myprog.c -lm
---
>ex12a_pth.px 10..1000000 1
np=1 nt=1 nx=10      it=30 rc=0 time=2.110936e-04 dmax=2.830711e-03
np=1 nt=1 nx=100     it=24 rc=0 time=2.509542e-04 dmax=2.869907e-05
np=1 nt=1 nx=1000    it=20 rc=0 time=6.354079e-04 dmax=2.871838e-07
np=1 nt=1 nx=10000   it=16 rc=0 time=4.451316e-03 dmax=1.578371e-08
np=1 nt=1 nx=100000  it=12 rc=0 time=3.547638e-02 dmax=9.606953e-07
np=1 nt=1 nx=1000000 it=7  rc=0 time=2.275494e-01 dmax=2.028484e-04
---
>ex12a_pth.px 100000 1..32
np=1 nt=1  nx=100000 it=12 rc=0 time=3.564769e-02 dmax=9.606953e-07
np=1 nt=2  nx=100000 it=12 rc=0 time=3.106295e-02 dmax=9.604200e-07
np=1 nt=4  nx=100000 it=12 rc=0 time=2.416283e-02 dmax=9.596533e-07
np=1 nt=8  nx=100000 it=12 rc=0 time=1.922178e-02 dmax=9.591225e-07
np=1 nt=16 nx=100000 it=12 rc=0 time=1.766124e-02 dmax=9.591722e-07
np=1 nt=32 nx=100000 it=12 rc=0 time=1.555847e-02 dmax=9.574738e-07
---
>mpirun -np 1..32 ex12a_pth.px 100000 1
np=1  nt=1 nx=100000 it=12 rc=0 time=5.037431e-02 dmax=9.606953e-07
np=2  nt=1 nx=100000 it=12 rc=0 time=3.466035e-02 dmax=9.604200e-07
np=4  nt=1 nx=100000 it=12 rc=0 time=1.898419e-02 dmax=9.596533e-07
np=8  nt=1 nx=100000 it=12 rc=0 time=1.825015e-02 dmax=9.591225e-07
np=16 nt=1 nx=100000 it=12 rc=0 time=1.199018e-02 dmax=9.591722e-07
np=32 nt=1 nx=100000 it=12 rc=0 time=7.261634e-03 dmax=9.574738e-07
---
>mpirun -np 1..32 ex12a_pth.px 100000 1..32
np=1  nt=1  nx=100000 it=12 rc=0 time=4.825157e-02 dmax=9.606953e-07
-
np=1  nt=2  nx=100000 it=12 rc=0 time=3.909390e-02 dmax=9.604200e-07
np=2  nt=1  nx=100000 it=12 rc=0 time=3.993990e-02 dmax=9.604200e-07
-
np=1  nt=4  nx=100000 it=12 rc=0 time=2.416283e-02 dmax=9.596533e-07
np=2  nt=2  nx=100000 it=12 rc=0 time=2.322900e-02 dmax=9.596533e-07
np=4  nt=1  nx=100000 it=12 rc=0 time=1.898419e-02 dmax=9.596533e-07
-
np=1  nt=8  nx=100000 it=12 rc=0 time=1.922178e-02 dmax=9.591225e-07
np=2  nt=4  nx=100000 it=12 rc=0 time=2.461819e-02 dmax=9.591225e-07
np=4  nt=2  nx=100000 it=12 rc=0 time=1.380887e-02 dmax=9.591225e-07
np=8  nt=1  nx=100000 it=12 rc=0 time=1.825015e-02 dmax=9.591225e-07
-
np=1  nt=16 nx=100000 it=12 rc=0 time=1.766124e-02 dmax=9.591722e-07
np=2  nt=8  nx=100000 it=12 rc=0 time=2.550384e-02 dmax=9.591722e-07
np=4  nt=4  nx=100000 it=12 rc=0 time=1.304034e-02 dmax=9.591722e-07
np=8  nt=2  nx=100000 it=12 rc=0 time=1.190888e-02 dmax=9.591722e-07
np=16 nt=1  nx=100000 it=12 rc=0 time=1.199018e-02 dmax=9.591722e-07

-
np=1  nt=32 nx=100000 it=12 rc=0 time=1.555847e-02 dmax=9.574738e-07
np=2  nt=16 nx=100000 it=12 rc=0 time=3.290873e-02 dmax=9.574738e-07
np=4  nt=8  nx=100000 it=12 rc=0 time=1.009655e-02 dmax=9.574738e-07
np=8  nt=4  nx=100000 it=12 rc=0 time=7.215660e-03 dmax=9.574738e-07
np=16 nt=2  nx=100000 it=12 rc=0 time=7.986512e-03 dmax=9.574738e-07
np=32 nt=1  nx=100000 it=12 rc=0 time=7.261634e-03 dmax=9.574738e-07
---