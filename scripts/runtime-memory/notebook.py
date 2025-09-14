#!/usr/bin/env python
# coding: utf-8

# In[2]:

from types import SimpleNamespace
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import numpy as np
import pandas as pd
import json, os, pysam, math, subprocess

# paths
path_workdir = Path("/home/b05b01002/HDD/project_nanoprep_re")
path_output = Path("./output")
os.makedirs(path_output, exist_ok=True)


# Find mprof.txt

# In[4]:


results = subprocess.run(
    f"find {path_workdir}/outputs -name mprof.txt".split(" "),
    capture_output=True, text=True, check=True
).stdout


# In[15]:


path_mprof = results.split("\n")
runtime_memory = pd.DataFrame()
for fname in path_mprof:
    if len(fname) == 0:
        continue
    fields = fname.split("/")
    software, beta_or_backend = fields[6], fields[8]
    data = pd.read_csv(fname, sep=" ")
    data = data.loc[data["CMDLINE"] == "MEM"]
    data = data.dropna(axis=1)
    data.columns = ["tag", "memory", "second"]
    data = data.loc[:, ["memory", "second"]]
    # down-sample to memory per second
    data["second"] = data["second"].apply(round)
    data.drop_duplicates(subset="second", keep="first", inplace=True, ignore_index=True)
    data["second"] = data["second"] - data.loc[0, "second"]
    data = data.loc[data["second"] % 5 == 0]
    data["memory"] = data["memory"] * 1.048576 / 1000 # MiB to GB
    data["software"] = software
    data["beta_or_backend"] = beta_or_backend if software == "Pychopper" else software
    runtime_memory = pd.concat([runtime_memory, data])

# In[16]:

fig, axes = plt.subplots(2, 1, figsize=(8, 10))
sns.lineplot(
    data=runtime_memory,
    x="second",
    y="memory",
    hue="beta_or_backend",
    errorbar="sd",
    estimator="mean",
    ax=axes[0]
)
sns.lineplot(
    data=runtime_memory[runtime_memory["second"] < 6000],
    x="second",
    y="memory",
    hue="beta_or_backend",
    errorbar="sd",
    estimator="mean",
    ax=axes[1]
)

fig.savefig("output/runtime-memory.svg")

