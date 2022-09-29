Installation instructions

This project is set up to be extendable to various SAT backends.
Only CaDiCaL is currently supported.

0) Go to the repository directory via

cd scm_sat

1) Link CaDiCaL using cmake via

cmake -DCMAKE_PREFIX_PATH="/path/to/cadical/install/dir/" .

e.g.

cmake -DCMAKE_PREFIX_PATH="/opt/cadical/" .

2) build project via

make 

3) call satscm binary with *constant*, *timeout* and *quiet* as arguments

./satscm *constant* *timeout* *quiet*

e.g.

./satscm 699839 300 1

to start the program for constant 699839 with 300 seconds timeout and without any debug messages.
Depending on your machine, the program should run for a few seconds (maybe up to a minute if your computer is slow?).
Once a solution is found, sat_scm prints a description of the resulting SCM circuit on the console.
