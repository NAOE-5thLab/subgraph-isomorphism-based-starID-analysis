gfortran -c detector/subgraph.f90
f2py --fcompiler=gnu95 -m subgraph_f -c --f90flags='-O3' detector/subgraph.f90
