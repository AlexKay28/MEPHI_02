IMM10:
---
>mpicc -o ex12a_mod.px -O2 ex12a_mod.c mycom.c mynet.c myprog.c -lm
---
>ex12a_mod.px 10..1000000
np=1 nx=10      it=30 rc=0 time=2.584606e-05 dmax=2.830711e-03
np=1 nx=100     it=24 rc=0 time=1.599714e-04 dmax=2.869907e-05
np=1 nx=1000    it=20 rc=0 time=1.094531e-03 dmax=2.871838e-07
np=1 nx=10000   it=16 rc=0 time=1.050541e-02 dmax=1.578371e-08
np=1 nx=100000  it=12 rc=0 time=4.193125e-02 dmax=9.606953e-07
np=1 nx=1000000 it=7  rc=0 time=2.244632e-01 dmax=2.028484e-04
---
>mpirun -np 1..32 ex12a_mod.px 100000
np=1  nx=100000 it=12 rc=0 time=4.951346e-02 dmax=9.606953e-07
np=2  nx=100000 it=12 rc=0 time=3.889497e-02 dmax=9.604200e-07
np=4  nx=100000 it=12 rc=0 time=2.767020e-02 dmax=9.596533e-07
np=8  nx=100000 it=12 rc=0 time=1.732946e-02 dmax=9.591225e-07
np=16 nx=100000 it=12 rc=0 time=1.238867e-02 dmax=9.591722e-07
np=32 nx=100000 it=12 rc=0 time=7.562608e-03 dmax=9.574738e-07
---
