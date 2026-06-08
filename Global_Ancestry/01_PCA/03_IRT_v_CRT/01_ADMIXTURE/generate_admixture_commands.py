import os
import glob

# User-defined parameters
input_dir = "/mnt/ceph/well0766/UMT_Redband/05_global_ancestry/01_STACKS/01_PCA/03_IRT_v_CRT/01_ADMIXTURE"
output_dir = "/mnt/ceph/well0766/UMT_Redband/05_global_ancestry/01_STACKS/01_PCA/03_IRT_v_CRT/01_ADMIXTURE/Unsupervised"
k_min = 2
k_max = 2
num_replicates = 40

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Find all .bed files
bed_files = glob.glob(os.path.join(input_dir, "*.bed"))

# Create the job file
job_file_path = os.path.join(output_dir, "admixture_jobs.sh")
with open(job_file_path, "w") as f:
    for bed_file in bed_files:
        base_name = os.path.splitext(bed_file)[0]
        prefix = os.path.basename(base_name)

        for k in range(k_min, k_max + 1):
            for rep in range(1, num_replicates + 1):
                seed = 1000 + (k * 100) + rep

                q_out = os.path.join(output_dir, f"{prefix}_K{k}_rep{rep}.Q")
                p_out = q_out.replace(".Q", ".P")
                cv_out = os.path.join(output_dir, f"{prefix}_K{k}_rep{rep}.cv")
                log_out = os.path.join(output_dir, f"{prefix}_K{k}_rep{rep}.log")

                cmd = (
                    f"cd {output_dir} && "  # Work from output directory to manage output files
                    f"admixture --cv --seed {seed} {bed_file} {k} > {log_out} 2>&1 && "
                    f"mv {prefix}.{k}.Q {q_out} && "
                    f"mv {prefix}.{k}.P {p_out} && "
                    f"mv log{prefix}.{k}.out {cv_out}"
                )
                f.write(cmd + "\n")

print(f"Job file written to: {job_file_path}")

