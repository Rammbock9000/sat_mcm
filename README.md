Installation instructions

This project is set up to be extendable to various SAT backends.
Only CaDiCaL is currently supported.

0) Go to the repository directory via

cd sat_scm

1) Link backends using cmake via

cmake -DCMAKE_PREFIX_PATH="/path/to/backend/install/dir/" .

CaDiCaL and Z3 are currently supported.

e.g.

cmake -DCMAKE_PREFIX_PATH="/opt/cadical/;/opt/z3/" .

2) build project via

make

3) call satscm binary with "constant", "backend name", "timeout" and "quiet" as arguments

./satscm "constant" "backend name" "timeout" "quiet"

e.g.

./satscm 14709 cadical 300 1

to start the program for constant 14709 with 300 seconds timeout and without any debug messages using cadical as backend solver.
Depending on your machine, the program should run for a few seconds (maybe up to a minute if your computer is slow?).
Once a solution is found, sat_scm prints a description of the resulting SCM circuit on the console.
