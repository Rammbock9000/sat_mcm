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

3) call satscm binary with "constant(s)", "backend name", "timeout", "threads" and "quiet" as arguments
-> specify multiple constants (for MCM instead of SCM) as a colon-separated list (e.g., 3:11:-9)
-> preprocessing steps for constants are
  1) performing right-shifts until it is an odd number
  2) taking the absolute value to avoid negative numbers
  3) removing duplicates

./satscm "constant(s)" "backend name" "timeout" "threads" "quiet"

e.g.

./satscm 14709 z3 300 2 1
to start the program for constant 14709 with 300 seconds timeout, 2 CPU threads allowed and without any debug messages using Z3 as backend solver.

./satscm 1:3:11:3:22 cadical 42 1 0
to start the program for constants 1, 3, 11, 3 and 22 (preprocessed to 3 and 11) with 42 seconds timeout, 1 CPU thread allowed and with debug messages using CaDiCaL as backend solver.

Note that CaDiCaL only supports single-thread solving.
Once a solution is found, sat_scm prints a description of the resulting SCM circuit on the console.
