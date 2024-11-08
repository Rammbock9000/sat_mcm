import subprocess


completed = set([
    "cmm",  # SKIP IT FOR NOW
    "pcmm", # SKIP IT FOR NOW
    "complex",
    "pcomplex",
    ("cmm", 4, 2, 2),
    ("cmm", 4, 2, 3),
    ("cmm", 4, 2, 4),
    ("cmm", 4, 2, 5),
    ("cmm", 4, 2, 6),
    ("cmm", 4, 2, 7),
    ("cmm", 4, 2, 8),
    ("cmm", 4, 2, 9),
    ("cmm", 4, 3, 2),
    ("cmm", 4, 3, 3),
    ("cmm", 4, 3, 4),
    ("cmm", 4, 3, 5),
    ("cmm", 4, 3, 6),
    ("cmm", 4, 3, 7),
    ("cmm", 4, 3, 8),
    ("cmm", 4, 3, 9),
    ("cmm", 4, 4, 2),
    ("cmm", 4, 4, 3),
    ("cmm", 4, 4, 4),
    ("cmm", 4, 4, 5),
    ("cmm", 4, 4, 6),
    ("cmm", 4, 4, 7),
    ("cmm", 4, 4, 8),
    ("cmm", 4, 4, 9),
    ("cmm", 4, 5, 2),
    ("cmm", 4, 5, 3),
    ("cmm", 4, 5, 4),
    ("cmm", 4, 5, 5),
    ("cmm", 4, 5, 6),
    ("cmm", 4, 5, 7),
    ("cmm", 4, 6, 2),
    ("cmm", 4, 6, 3),
    ("cmm", 4, 6, 4),
    ("cmm", 4, 6, 6),
    ("cmm", 4, 8, 3),
    ("cmm", 4, 9, 3),
    ("cmm", 6, 2, 2),
    ("cmm", 6, 2, 3),
    ("cmm", 6, 2, 4),
    ("cmm", 6, 3, 2),
    ("cmm", 6, 3, 3),
    ("cmm", 6, 4, 2),
    ("cmm", 6, 4, 3),
    ("cmm", 6, 5, 2),
    ("cmm", 6, 6, 2),
    ("cmm", 6, 7, 2),
    ("cmm", 6, 8, 2),
    ("cmm", 6, 9, 2),
    ("cmm", 8, 2, 2),
    ("cmm", 8, 2, 3),
    ("cmm", 8, 3, 2),
    ("cmm", 8, 4, 2),
    ("cmm", 8, 5, 2),
    ("cmm", 8, 6, 2),
    ("cmm", 8, 9, 2),
    ("conv", 2, 9),
])


def is_completed(experiment, W=None, M=None, N=None):
    if experiment in completed:
        return True
    #return (experiment, W, M, N) in completed
    return False


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
    filename = f"{base_dir}benchmark/linux-cluster/start_cmm-{experiment}_{additional_args.replace(' ', '_')}.sh"
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
    filename = f"{base_dir}{script_base}/start_cleanup-{experiment}_{get_additional_args(experiment, W, M, N).replace(' ', '_')}.sh"
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
    no_pipe_experiments = ["cmm", "conv", "complex"]
    pipe_experiments = ["pcmm", "pconv", "pcomplex"]
    experiments = no_pipe_experiments + pipe_experiments # pipe_experiments + no_pipe_experiments
    how_often = 1
    num_submitted_total = 0
    max_jobs = 498
    limit_reached = False
    DEBUGGING = False
    for experiment in experiments:
        num_submitted = 0
        if is_completed(experiment):
            continue
        Ws, Ms, Ns = get_loops(experiment)
        for W in Ws:
            for M in Ms:
                for N in Ns:
                    if limit_reached or is_completed(experiment, W, M, N):
                        continue
                    if not DEBUGGING:
                        clean_file = create_cleanup_slurm_script(experiment, W, M, N)
                        cmm_file = create_cmm_slurm_script(experiment, W, M, N)
                    job_id = None
                    for _ in range(how_often):
                        if num_submitted_total + num_submitted + 1 > max_jobs:
                            limit_reached = True
                            break
                        else:
                            if not DEBUGGING:
                                job_id = start_job_and_get_id(clean_file, dep=job_id)
                            num_submitted += 1
                        if num_submitted_total + num_submitted + 1 > max_jobs:
                            limit_reached = True
                            break
                        else:
                            if not DEBUGGING:
                                job_id = start_job_and_get_id(cmm_file, dep=job_id)
                            num_submitted += 1
        num_submitted_total += num_submitted
        print(f"Submitted {num_submitted} jobs for experiment '{experiment}'")
    print(f"Submitted {num_submitted_total} jobs in total!")


if __name__ == '__main__':
    main()