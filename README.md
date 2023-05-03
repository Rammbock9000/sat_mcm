Installation instructions

This project is set up to be extendable to various SAT backends.
CaDiCaL, Glucose-Syrup and Z3 are currently supported.

0) Go to the repository directory via

cd sat_mcm

1) Link backends using cmake via

cmake -DCMAKE_PREFIX_PATH="/path/to/backend/install/dir/" .

e.g.

cmake -DCMAKE_PREFIX_PATH="/opt/cadical/;/opt/z3/" .

2) build project via

make

3) call satmcm binary without arguments to see command line arguments.

Once a solution is found, sat_mcm prints a description of the resulting MCM circuit on the console.

4) (optional) install everything (to default location) via

make install

Step 4 likely makes sense if you want to use our project as a library. 

You can also specify another install dir in step 1) via

cmake -DCMAKE_PREFIX_PATH="/path/to/backend/install/dir/" -DCMAKE_INSTALL_PREFIX="/path/to/install/dir/" .

e.g.

cmake -DCMAKE_PREFIX_PATH="/opt/cadical/;/opt/z3/" -DCMAKE_INSTALL_PREFIX="/opt/mcm/" .

To reproduce all our experiments do the following:

1) get a working python 3 environment

2) go to the repository root directory (where e.g. the CMakeLists.txt and the binary are located)

3) call "python benchmark/scripts/start_benchmark.py" (this took approx. 1 month of compute time on our machines...)

4) call "python benchmark/scripts/create_vhdl_code.py" (this finishes within a few seconds)
