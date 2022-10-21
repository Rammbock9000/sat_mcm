Installation instructions

This project is set up to be extendable to various SAT backends.
CaDiCaL and Z3 are currently supported.

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
  i) performing right-shifts until it is an odd number
  ii) taking the absolute value to avoid negative numbers
  iii) removing duplicates

./satscm "constant(s)" "backend name" "timeout" "threads" "quiet" "optimize full adders"

e.g.

./satscm 1:3:11:3:22 cadical 42 1 0 1
to start the program for constants 1, 3, 11, 3 and 22 (preprocessed to 3 and 11) with 42 seconds timeout, 1 CPU thread allowed and with debug messages using CaDiCaL as backend solver and optimizing the full adder count for the optimal number of adders.

Note that CaDiCaL only supports single-thread solving.
Once a solution is found, sat_scm prints a description of the resulting SCM circuit on the console.
