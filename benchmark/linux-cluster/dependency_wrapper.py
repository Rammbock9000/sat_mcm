import subprocess


def create_cmm_slurm_script(experiment):
    time_limit = f"6-00:00:00"
    num_tasks = 48
    time_limit = f"0-00:01:00"
    num_tasks = 2
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
        f.write(f"#SBATCH --output={base_dir}benchmark/linux-cluster/cmm-{experiment}.out\n")
        if experiment == "cmm":
            f.write(f"cd {base_dir};python3 {script_base}/cmm_benchmark.py {num_tasks-1}\n")
        elif experiment == "pcmm":
            f.write(f"cd {base_dir};python3 {script_base}/pcmm_benchmark.py {num_tasks-1}\n")
        elif experiment == "complex":
            f.write(f"cd {base_dir};python3 {script_base}/complex_mult_benchmark.py {num_tasks-1}\n")
        elif experiment == "pcomplex":
            f.write(f"cd {base_dir};python3 {script_base}/pcomplex_mult_benchmark.py {num_tasks-1}\n")
        elif experiment == "conv":
            f.write(f"cd {base_dir};python3 {script_base}/conv_core_benchmark.py {num_tasks-1}\n")
        elif experiment == "pconv":
            f.write(f"cd {base_dir};python3 {script_base}/pconv_core_benchmark.py {num_tasks-1}\n")
    return filename


def create_cleanup_slurm_script(experiment):
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
        f.write(f"#SBATCH --time=0-00:10:00\n")
        f.write(f"#SBATCH --mail-user=nfiege@uni-kassel.de\n")
        f.write(f"#SBATCH --mail-type=ALL\n")
        f.write(f"#SBATCH --output={base_dir}benchmark/linux-cluster/clean-{experiment}.out\n")
        f.write(f"cd ${base_dir};python3 {script_base}/cleanup_intermediate.py {result_dir} {experiment}\n")
    return filename


def start_job_and_get_id(job, dep=None):
    if dep is None:
        command = ["sbatch", f"{job}"]
    else:
        command = ["sbatch", f"--dependency=afterany:{dep}", job]
    print(f"PYTHON: start executing command: {command}")
    output = subprocess.run(command, stdout=subprocess.PIPE).stdout.decode('utf-8').rstrip()
    print(f"PYTHON: finished")
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


def main():
    experiments = ["cmm", "pcmm", "conv", "pconv", "complex", "pcomplex"]
    experiments = ["cmm"]
    for experiment in experiments:
        cmm_file = create_cmm_slurm_script(experiment)
        clean_file = create_cleanup_slurm_script(experiment)
        job_id_1 = start_job_and_get_id(cmm_file)
        job_id_2 = start_job_and_get_id(clean_file, dep=job_id_1)
        job_id_3 = start_job_and_get_id(cmm_file, dep=job_id_2)
        #job_id_4 = start_job_and_get_id(clean_file, dep=job_id_3)
        #job_id_5 = start_job_and_get_id(cmm_file, dep=job_id_4)
        #job_id_6 = start_job_and_get_id(clean_file, dep=job_id_5)
        #job_id_7 = start_job_and_get_id(cmm_file, dep=job_id_6)
        #_ = start_job_and_get_id(clean_file, dep=job_id_7)


if __name__ == '__main__':
    main()