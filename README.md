## General Information

Contents of the main branch are published under:

N. Fiege, M. Kumm and P. Zipf, "Bit-Level Optimized Constant Multiplication Using Boolean Satisfiability," in IEEE Transactions on Circuits and Systems I: Regular Papers, vol. 71, no. 1, pp. 249-261, Jan. 2024, doi: 10.1109/TCSI.2023.3327814.

Please cite our work if you found it helpful:

    @ARTICLE{10320216,
        author={Fiege, Nicolai and Kumm, Martin and Zipf, Peter},
        journal={IEEE Transactions on Circuits and Systems I: Regular Papers},
        title={Bit-Level Optimized Constant Multiplication Using Boolean Satisfiability},
        year={2024},
        volume={71},
        number={1},
        pages={249-261},
        doi={10.1109/TCSI.2023.3327814}
    }
    
Contents of the CMM branch are currently under review for publication... Please cite the above paper until then.

## Installation instructions

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

## Reproduction of SCM/MCM results

To reproduce all our experiments do the following:

1) get a working python 3 environment

2) go to the repository root directory (where e.g. the CMakeLists.txt and the binary are located)

3) call "python benchmark/scripts/start_benchmark.py" to start the experiments (this took approx. 1 month of compute time on our machines...)

4) call "python benchmark/scripts/create_vhdl_code.py" to generate vhdl code (this finishes within a few seconds)

5) call "python benchmark/scripts/analyze.py" to analyze experimental results (this might take a few minutes)

6) call "python benchmark/scripts/synthesize_vivado.py" to generate synthesis results with vivado (this might take an hour or two)

Now, you should find pre-synthesis results summarized in "benchmark/results/results.csv" and post-place-and-route results (timing, area, power consumption) in the respective vivado project directories at location "benchmark/vhdl/*experiment_type*/*setting*/*vivado_project*/"
