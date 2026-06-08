import os
import glob
import numpy as np

# User parameters
input_dir = "."  # Change this if needed
output_avg_q = os.path.join(input_dir, "averaged_admixture_Q.txt")
output_std_q = os.path.join(input_dir, "admixture_stddev_Q.txt")
output_flags = os.path.join(input_dir, "low_minor_components.txt")
output_avg_p = os.path.join(input_dir, "averaged_admixture_P.txt")

# Thresholds for flagging minor ancestry
low_ancestry_threshold = 0.01
instability_threshold = 0.005

# Step 1: Locate all .Q and .P files (paired by prefix)
q_files = sorted(glob.glob(os.path.join(input_dir, "*.Q")))
p_files = sorted(glob.glob(os.path.join(input_dir, "*.P")))

if not q_files:
    raise ValueError(f"No .Q files found in {input_dir}")

# Step 2: Load all Q files
q_matrices = [np.loadtxt(qf) for qf in q_files]
q_shapes = [m.shape for m in q_matrices]
if len(set(q_shapes)) != 1:
    raise ValueError(f"Not all Q files have the same dimensions: {set(q_shapes)}")

stacked_q = np.stack(q_matrices, axis=0)
average_q = np.mean(stacked_q, axis=0)
std_dev_q = np.std(stacked_q, axis=0)

# Step 3: Save averaged Q and std dev Q
np.savetxt(output_avg_q, average_q, fmt="%.6f", delimiter=' ')
np.savetxt(output_std_q, std_dev_q, fmt="%.6f", delimiter=' ')
print(f"Averaged Q saved to: {output_avg_q}")
print(f"Standard deviation Q saved to: {output_std_q}")

# Step 4: Flag low minor ancestry components
with open(output_flags, "w") as out:
    out.write("Individual\tComponent(K)\tMean\tStdDev\n")
    for i in range(average_q.shape[0]):  # individuals
        for j in range(average_q.shape[1]):  # clusters
            mean = average_q[i, j]
            std = std_dev_q[i, j]
            if mean < low_ancestry_threshold and std > instability_threshold:
                out.write(f"Ind_{i+1}\tK{j+1}\t{mean:.6f}\t{std:.6f}\n")

print(f"Flagged low minor components saved to: {output_flags}")

# Step 5 (Optional): Average .P files if available
if p_files:
    p_matrices = [np.loadtxt(pf) for pf in p_files]
    p_shapes = [m.shape for m in p_matrices]
    if len(set(p_shapes)) != 1:
        raise ValueError(f"Not all P files have the same dimensions: {set(p_shapes)}")

    stacked_p = np.stack(p_matrices, axis=0)
    average_p = np.mean(stacked_p, axis=0)
    np.savetxt(output_avg_p, average_p, fmt="%.6f", delimiter=' ')
    print(f"Averaged P saved to: {output_avg_p}")
else:
    print("No .P files found — skipping P averaging.")
