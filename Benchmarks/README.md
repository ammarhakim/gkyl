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

