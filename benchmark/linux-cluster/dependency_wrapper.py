import subprocess


def get_additional_args(experiment, W, M, N):
    if "complex" in experiment:
        return f"{W}"
    if "conv" in experiment:
        return f"{W} {N}"
    # "conv"/"pconv"
    return f"{W} {M} {N}"


def create_cmm_slurm_script(experiment, W, M, N):
    additional_args = get_additional_args(experiment, W, M, N)
    time_limit = f"6-00:00:00"
    num_tasks = 48
    mem = 5000*num_tasks
    base_dir="/home/groups/fb16-digi/projects/constmult/sat_mcm/"
    script_base="benchmark/scripts"
    filename = f"{base_dir}benchmark/linux-cluster/start_cmm-{experiment}.sh"
    with open(filename, "w") as f:
        f.write(f"#!/bin/bash\n")
        f.write(f"#SBATCH --nodes=1\n")
        f.write(f"#SBATCH --tasks-per-node={num_tasks}\n")
        f.write(f"#SBATCH -p pub23\n")
        f.write(f"#SBATCH --mem={mem}\n")
        f.write(f"#SBATCH --time={time_limit}\n")
        f.write(f"#SBATCH --mail-user=nfiege@uni-kassel.de\n")
        f.write(f"#SBATCH --mail-type=ALL\n")
        f.write(f"#SBATCH --output={base_dir}benchmark/linux-cluster/cmm-{experiment}_{additional_args.replace(' ', '_')}.out\n")
        if experiment == "cmm":
            f.write(f"cd {base_dir};python3 {script_base}/cmm_benchmark.py {num_tasks-1} {additional_args}\n")
        elif experiment == "pcmm":
            f.write(f"cd {base_dir};python3 {script_base}/pcmm_benchmark.py {num_tasks-1} {additional_args}\n")
        elif experiment == "complex":
            f.write(f"cd {base_dir};python3 {script_base}/complex_mult_benchmark.py {num_tasks-1} {additional_args}\n")
        elif experiment == "pcomplex":
            f.write(f"cd {base_dir};python3 {script_base}/pcomplex_mult_benchmark.py {num_tasks-1} {additional_args}\n")
        elif experiment == "conv":
            f.write(f"cd {base_dir};python3 {script_base}/conv_core_benchmark.py {num_tasks-1} {additional_args}\n")
        elif experiment == "pconv":
            f.write(f"cd {base_dir};python3 {script_base}/pconv_core_benchmark.py {num_tasks-1} {additional_args}\n")
        else:
            raise Exception(f"invalid experiment requested: '{experiment}'")
    return filename


def create_cleanup_slurm_script(experiment, W, M, N):
    base_dir="/home/groups/fb16-digi/projects/constmult/sat_mcm/"
    script_base="benchmark/linux-cluster"
    filename = f"{base_dir}{script_base}/start_cleanup-{experiment}.sh"
    result_dir=f"{base_dir}benchmark/results/cmm_rnd"
    with open(filename, "w") as f:
        f.write(f"#!/bin/bash\n")
        f.write(f"#SBATCH --nodes=1\n")
        f.write(f"#SBATCH --tasks-per-node=1\n")
        f.write(f"#SBATCH -p pub23\n")
        f.write(f"#SBATCH --mem=4000\n")
        f.write(f"#SBATCH --time=0-00:30:00\n")
        f.write(f"#SBATCH --output={base_dir}benchmark/linux-cluster/clean-{experiment}_{get_additional_args(experiment, W, M, N).replace(' ', '_')}.out\n")
        f.write(f"cd {base_dir};python3 {script_base}/cleanup_intermediate.py {result_dir} {experiment} {W} {M} {N}\n")
    return filename


def start_job_and_get_id(job, dep=None):
    if dep is None:
        command = ["sbatch", f"{job}"]
    else:
        command = ["sbatch", f"--dependency=afterany:{dep}", job]
    output = subprocess.run(command, stdout=subprocess.PIPE).stdout.decode('utf-8').rstrip()
    elems = output.split()
    if len(elems) < 4:
        raise Exception(f"failed to queue batch jobs (code 0 - '{output}')")
    if elems[0] != "Submitted":
        raise Exception(f"failed to queue batch jobs (code 1 - '{elems[0]}')")
    if elems[1] != "batch":
        raise Exception(f"failed to queue batch jobs (code 2 - '{elems[1]}')")
    if elems[2] != "job":
        raise Exception(f"failed to queue batch jobs (code 3 - '{elems[2]}')")
    return elems[3]


def get_loops(experiment):
    # Ws, Ms, Ns
    if "complex" in experiment:
        return [2,3,4,5,6,7,8,9], [-1], [-1]
    if "conv" in experiment:
        return [2,4,6,8], [-1], [2,3,4,5,6,7,8,9]
    # "conv"/"pconv"
    return [4,6,8], [2,3,4,5,6,7,8,9], [2,3,4,5,6,7,8,9]


def main():
    experiments = ["cmm", "pcmm", "conv", "pconv", "complex", "pcomplex"]
    how_often = 1
    for experiment in experiments:
        Ws, Ms, Ns = get_loops(experiment)
        for W in Ws:
            for M in Ms:
                for N in Ns:
                    cmm_file = create_cmm_slurm_script(experiment, W, M, N)
                    clean_file = create_cleanup_slurm_script(experiment, W, M, N)
                    job_id = start_job_and_get_id(clean_file)
                    for _ in range(how_often):
                        job_id = start_job_and_get_id(cmm_file, dep=job_id)
                        job_id = start_job_and_get_id(clean_file, dep=job_id)


if __name__ == '__main__':
    main()