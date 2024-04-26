import subprocess


def start_job_and_get_id(job, dep=None):
    if dep is None:
        command = ["sbatch", f"{job}"]
    else:
        command = ["sbatch", f"--dependency=afterany:{dep} {job}"]
    output = subprocess.run(command, stdout=subprocess.PIPE).stdout.decode('utf-8').rstrip()
    elems = output.split()
    if len(elems) < 5:
        raise Exception("failed to queue batch jobs (code 0)")
    if elems[0] != "successfully":
        raise Exception("failed to queue batch jobs (code 1)")
    if elems[1] != "submitted":
        raise Exception("failed to queue batch jobs (code 2)")
    if elems[2] != "batch":
        raise Exception("failed to queue batch jobs (code 3)")
    if elems[3] != "job":
        raise Exception("failed to queue batch jobs (code 4)")
    return elems[4]


def main():
    experiments = ["cmm", "pcmm", "conv", "pconv", "complex", "pcomplex"]
    for experiment in experiments:
        job_id_1 = start_job_and_get_id(f"start_cmm.sh {experiment}")
        job_id_2 = start_job_and_get_id(f"start_cleanup.sh {experiment}", dep=job_id_1)
        job_id_3 = start_job_and_get_id(f"start_cmm.sh {experiment}", dep=job_id_2)
        job_id_4 = start_job_and_get_id(f"start_cleanup.sh {experiment}", dep=job_id_3)
        job_id_5 = start_job_and_get_id(f"start_cmm.sh {experiment}", dep=job_id_4)
        job_id_6 = start_job_and_get_id(f"start_cleanup.sh {experiment}", dep=job_id_5)
        job_id_7 = start_job_and_get_id(f"start_cmm.sh {experiment}", dep=job_id_6)
        _ = start_job_and_get_id(f"start_cleanup.sh {experiment}", dep=job_id_7)


if __name__ == '__main__':
    main()