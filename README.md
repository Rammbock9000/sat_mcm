Installation instructions

This project is set up to be extendable to various SAT backends.
CaDiCaL, Glucose-Syrup and Z3 are currently supported.

0) Go to the repository directory via

cd sat_scm

1) Link backends using cmake via

cmake -DCMAKE_PREFIX_PATH="/path/to/backend/install/dir/" .

e.g.

cmake -DCMAKE_PREFIX_PATH="/opt/cadical/;/opt/z3/" .

2) build project via

make

3) call satscm binary without arguments to see command line arguments.

Once a solution is found, sat_scm prints a description of the resulting SCM circuit on the console.
