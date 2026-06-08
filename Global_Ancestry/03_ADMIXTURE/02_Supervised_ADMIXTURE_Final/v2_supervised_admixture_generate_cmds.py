import os
import glob

# ============================
# User-defined parameters
# ============================
input_dir = "/mnt/ceph/well0766/UMT_Redband/05_global_ancestry/01_STACKS/03_ADMIXTURE/01_input_files/03_Generate_Plink_Files"
output_dir = "/mnt/ceph/well0766/UMT_Redband/05_global_ancestry/01_STACKS/03_ADMIXTURE/02_Supervised_ADMIXTURE_Final/Results"
k_min = 4
k_max = 4
num_replicates = 100

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Find all .bed files
bed_files = glob.glob(os.path.join(input_dir, "*.bed"))

# ============================
# Create the job file
# ============================
job_file_path = os.path.join(output_dir, "admixture_supervised_jobs.sh")
with open(job_file_path, "w") as f:
    for bed_file in bed_files:
        prefix = os.path.splitext(os.path.basename(bed_file))[0]

        for k in range(k_min, k_max + 1):
            for rep in range(1, num_replicates + 1):
                seed = 1000 + (k * 100) + rep

                # Temporary working directory per replicate
                tmp_dir = os.path.join(output_dir, f"{prefix}_K{k}_rep{rep}_tmp")

                # Output files
                q_out = os.path.join(output_dir, f"{prefix}_K{k}_rep{rep}.Q")
                p_out = os.path.join(output_dir, f"{prefix}_K{k}_rep{rep}.P")
                log_out = os.path.join(output_dir, f"{prefix}_K{k}_rep{rep}.log")

                cmd = (
                    f"mkdir -p {tmp_dir} && "
                    f"cd {tmp_dir} && "
                    f"admixture --supervised --seed {seed} {bed_file} {k} > {log_out} 2>&1 && "
                    f"mv {prefix}.{k}.Q {q_out} && "
                    f"mv {prefix}.{k}.P {p_out} && "
                    f"cd {output_dir} && "
                    f"rm -rf {tmp_dir}"
                )
                
                f.write(cmd + "\n")

print(f"✅ Supervised job file written to: {job_file_path}")
