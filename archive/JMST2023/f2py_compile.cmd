gfortran -c matching/subgraph.f90
f2py --fcompiler=gnu95 --compiler=mingw32 -m subgraph_f -c --f90flags='-O3' matching/subgraph.f90
