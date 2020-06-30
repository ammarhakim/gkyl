# Tests for benchmarking

This directory contains tests used to benchmark the code on various
platforms. Some of these tests take significant time to run on modest
(or even large) number of cores. They are mainly designed to
understand performace characteristics of the machines which Gkeyll
runs on.

In each directory there are sub-directories with log-file from
specific machines. If you run the benchmark on a new machine please
add details for that machine and the log-file to that directory.

# vm-weibel-4d: Vlasov-LBO-Maxwell System (Weibel instability)

4D Weibel instability problem. Solves the VM system with LBO
collisions. Ions are not evolved. See ApJ Letters (2019) for details
of the physics studied.

# oblique-mode-4d: Vlasov-Maxwell System (Oblique, Hybrid Two-Stream-Weibel, mode)

4D oblique, hybrid two-stream-weibel, mode simulation. Solves the VM system with
no collision operator. Ions are not evolved. The no-app input file in the portal-V100
folder can be run on a GPU.

Performance benchmarks for 1000 time-steps
Nvidia V100 GPU - 884.368 seconds
Stampede 2 Intel Skylake (Intel Xeon Platinum 8160, 48 cores per node) - 1888.513 seconds
Frontera Intel Cascade Lake (Intel Xeon Platinum 8280, 56 cores per node) - 1697.106 seconds

