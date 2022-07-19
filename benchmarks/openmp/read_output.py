from glob import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

results = []
for res_file in glob("results/*/omp.csv"):
    df = pd.read_csv(res_file)
    df["filename"] = res_file
    results.append(df)

df = pd.concat(results, axis="index")
df["commit_hash"] = df["filename"].str.extract("/(.*)/")
df["short_commit_hash"] = df["commit_hash"].str.slice(0, 6)

df.columns = df.columns.str.strip()
ref_time = df.loc[df["n_threads"] == 1]["elapsed_time"].loc[0]

df["speedup"] = ref_time/df["elapsed_time"]
df["efficiency"] = df["speedup"]/df["n_threads"]

plt.figure()
for (commit, group_df) in df.groupby("short_commit_hash"):
    plt.plot(group_df["n_threads"], group_df["speedup"], label=commit)
n_threads_range = np.arange(df.n_threads.min(), df.n_threads.max()+1)
plt.plot(n_threads_range, n_threads_range, linestyle="--", color="grey")
plt.xlabel("Number of threads")
plt.ylabel("Speedup")
plt.title("OpenMP parallelization")
plt.legend()
plt.show()
