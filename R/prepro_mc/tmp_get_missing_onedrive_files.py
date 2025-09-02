from glob import glob

path = "/Users/johanna/Uni/masterarbeit/output/mof3d_submontane/Results"
#path = "/Users/johanna/Uni/masterarbeit/output/lowland/Results"

missing_files = []

req_files_prefix = ["light", "voxel", "trees", "shoots", "mortality"]

n_files = 0

for prefix in req_files_prefix:
    files = glob(f"{path}/{prefix}*replicate_0_time_step*.txt", recursive=True)

    n_files += len(files)

    if len(files) == 0:
        missing_files += [f"{prefix}_replicate_0_time_step_{i}.txt" for i in range(1, 201)]
        print(f"All files missing for prefix: {prefix}")
    if len(files) < 201:
        existing_time_steps = [int(f.split("_")[-1].split(".")[0]) for f in files]
        new_missing = [f"{prefix}_replicate_0_time_step_{i}.txt"
                          for i in range(1, 201) if i not in existing_time_steps]
        missing_files += new_missing
        print(f"{len(new_missing)} files missing for prefix: {prefix}")
    else:
        print("All files found for prefix:", prefix)

for f in missing_files:
    print(f)

# ----- Second check -----

# Count lines in each file and print warning if it is 0
for prefix in req_files_prefix:
    files = glob(f"{path}/{prefix}*replicate_0_time_step*.txt", recursive=True)
    for file in files:
        with open(file, 'r') as f:
            lines = f.readlines()
            if len(lines) < 5:
                print(f"Warning: {file} is empty.")
