#!/usr/bin/env python
# coding: utf-8

# In[2]:


from types import SimpleNamespace
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import numpy as np
import pandas as pd
import json, os, pysam, math

# paths
path_workdir = Path("/home/b05b01002/HDD/project_nanoprep_re")
path_output = Path("./output")
os.makedirs(path_output, exist_ok=True)


# In[3]:


# NanoPrePro/Pychopper full-length read file name template
path_fastq = lambda wildcards: \
    path_workdir / f"outputs/{wildcards.pretty_name}/mprof_{wildcards.software}/{wildcards.beta_or_backend}/{wildcards.name}_{wildcards.accuracy}/full-length.fq"

# wildcards
names_accuracies = {
    "egr-109-bio1": ["sup", "hac", "fast"],
    "egr-109-bio2": ["sup", "hac", "fast"],
    "lch-109-bio1": ["sup", "hac", "fast"],
    "lch-109-bio2": ["sup", "hac", "fast"],
    "ptr-109-bio1": ["sup", "hac", "fast"],
    "ptr-109-bio2": ["sup", "hac", "fast"],
    "ptr-111-bio1": ["sup", "hac", "fast"],
    "ont-10x-human": ["sup", "hac", "fast"],
    "ont-visium-mouse": ["sup", "hac", "fast"],
    "mouse-retina-subset1": ["pre-called"],
    "mouse-retina-subset2": ["pre-called"]
}
accuracies = [
    "sup",
    "hac",
    "fast",
    "pre-called"
]
betas = [f"beta{i}" for i in ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0", "1.5", "2.0"]]
backends = ["phmm", "edlib"]
softwares = ["pychopper", "nanoprep"]
pretty_names = {
    "pychopper": "Pychopper",
    "nanoprep": "NanoPreP"
}


# Pychopper/NanoPreP processed read length

# In[4]:


data = pd.DataFrame()
for name, accuracies in names_accuracies.items():
    for accuracy in accuracies:
        for software in softwares:
            if (accuracy == "pre-called") & (software == "pychopper"):
                continue
            beta_or_backends = {
                "nanoprep": betas,
                "pychopper": backends
            }[software]
            clean_name = {
                "pychopper": lambda x: x.split("|")[-1],
                "nanoprep": lambda x: x
            }[software]
            for beta_or_backend in beta_or_backends:
                wildcards = SimpleNamespace(
                    pretty_name=pretty_names[software],
                    software=software,
                    beta_or_backend=beta_or_backend,
                    name=name,
                    accuracy=accuracy
                )
                fname = path_fastq(wildcards)
                with pysam.FastxFile(fname) as handle:
                    name_len = [[clean_name(i.name), len(i.sequence)] for i in handle]
                data_tmp = pd.DataFrame(name_len)
                data_tmp.columns = ["read_name", "length"]
                data_tmp["software"] = software
                data_tmp["beta_or_backend"] = beta_or_backend
                data_tmp["name"] = name
                data_tmp["accuracy"] = accuracy if not accuracy == "pre-called" else "sup"
                data = pd.concat([data, data_tmp])


# In[ ]:


data.to_csv("processed_length.csv", index=False)


# Raw read length

# In[ ]:


# raw read file name template
path_fastq = lambda wildcards: \
    path_workdir / f"outputs/Basecalling/aggregate_fastq/{wildcards.name}_{wildcards.accuracy}.fq"

# wildcards
names_accuracies = {
    "egr-109-bio1": ["sup", "hac", "fast"],
    "egr-109-bio2": ["sup", "hac", "fast"],
    "lch-109-bio1": ["sup", "hac", "fast"],
    "lch-109-bio2": ["sup", "hac", "fast"],
    "ptr-109-bio1": ["sup", "hac", "fast"],
    "ptr-109-bio2": ["sup", "hac", "fast"],
    "ptr-111-bio1": ["sup", "hac", "fast"],
    "ont-10x-human": ["sup", "hac", "fast"],
    "ont-visium-mouse": ["sup", "hac", "fast"],
}
precalled_fastq = {
    "mouse-retina-subset1": path_workdir / "rawdata/rui-chen-lab/Ms_bulk_subset1.fastq.gz",
    "mouse-retina-subset2": path_workdir / "rawdata/rui-chen-lab/Ms_bulk_subset2.fastq.gz"
}


# Basecalled fastq

# In[ ]:


data_raw = pd.DataFrame()
for name, accuracies in names_accuracies.items():
    for accuracy in accuracies:
        wildcards = SimpleNamespace(name=name, accuracy=accuracy)
        fname = path_fastq(wildcards)
        with pysam.FastxFile(fname) as handle:
            name_len = [[i.name, len(i.sequence)] for i in handle]
        data_tmp = pd.DataFrame(name_len)
        data_tmp.columns = ["read_name", "length"]
        data_tmp["software"] = "raw"
        data_tmp["beta_or_backend"] = "raw"
        data_tmp["name"] = name
        data_tmp["accuracy"] = accuracy
        data_raw = pd.concat([data_raw, data_tmp])


# Pre-called fastq

# In[ ]:


for name, fname in precalled_fastq.items():
    with pysam.FastxFile(fname) as handle:
        name_len = [[i.name, len(i.sequence)] for i in handle]
    data_tmp = pd.DataFrame(name_len)
    data_tmp.columns = ["read_name", "length"]
    data_tmp["software"] = "raw"
    data_tmp["beta_or_backend"] = "raw"
    data_tmp["name"] = name
    data_tmp["accuracy"] = "sup"
    data_raw = pd.concat([data_raw, data_tmp])


# In[ ]:


data_raw.to_csv("raw_length.csv", index=False)


# Combine processed and raw read length

# In[ ]:


data = pd.concat([
    pd.read_csv("processed_length.csv"),
    pd.read_csv("raw_length.csv")
])


# In[ ]:


n_samples = 10000
beta_or_backend_to_plot = {
    "nanoprep": betas,
    "pychopper": backends
}
names_accuracies = {
    "egr-109-bio1": ["sup", "hac", "fast"],
    "egr-109-bio2": ["sup", "hac", "fast"],
    "lch-109-bio1": ["sup", "hac", "fast"],
    "lch-109-bio2": ["sup", "hac", "fast"],
    "ptr-109-bio1": ["sup", "hac", "fast"],
    "ptr-109-bio2": ["sup", "hac", "fast"],
    "ptr-111-bio1": ["sup", "hac", "fast"],
    "ont-10x-human": ["sup", "hac", "fast"],
    "ont-visium-mouse": ["sup", "hac", "fast"],
    "mouse-retina-subset1": ["sup"],
    "mouse-retina-subset2": ["sup"]
}
# 
data_plot = pd.DataFrame()
for name, accuracies in names_accuracies.items():
    print(name)
    for accuracy in accuracies:
        print(accuracy, end=" ")
        is_name = data["name"] == name
        is_accuracy = data["accuracy"] == accuracy
        is_raw = data["beta_or_backend"] == "raw"
        data_x = data.loc[is_name & is_accuracy & is_raw]
        for software in softwares:
            if (software == "pychopper") & (name in ["mouse-retina-subset1", "mouse-retina-subset2"]):
                continue
            print(software, end=" ")
            for beta_or_backend in beta_or_backend_to_plot[software]:
                is_beta_or_backend = data["beta_or_backend"] == beta_or_backend
                data_y = data.loc[is_name & is_accuracy & is_beta_or_backend]
                data_y = data_y.sample(n=round(n_samples / len(beta_or_backend_to_plot[software])), replace=False)
                data_xy = pd.merge(data_x, data_y, how="right", on=["read_name"])
                data_plot = pd.concat([data_plot, data_xy])
                print(beta_or_backend, end=" ")
        print()


# In[ ]:


data_plot.to_csv("downsampled_plot_data.csv", index=False)


# In[ ]:


palette = {}
palette["phmm"], palette["edlib"] = sns.color_palette(["#8a6d53", "#ab7849"], 2)
for i, j in zip(betas, sns.color_palette("viridis_r", len(betas))):
    palette[i] = j

col_wrap = 4
col = "name_x"
col_values = data_plot[col].unique()
x="length_x"
y="length_y"
hue = "beta_or_backend_y"
style = "accuracy_x"
markers = {"sup": "o", "hac": "^", "fast": "s"}
subplot_aspect = 1
subplot_height = 4
subplot_width = subplot_height * subplot_aspect
fig_rows = math.ceil(len(col_values) / col_wrap)
fig_cols = col_wrap

# 
fig, axes = plt.subplots(fig_rows, col_wrap, figsize=(subplot_width * fig_cols, subplot_height * fig_rows))
plt.subplots_adjust(wspace=0.1, hspace=0.1)
for name, ax in zip(col_values, axes.flatten()):
    sns.scatterplot(
        data_plot[data_plot["name_x"] == name],
        x=x,
        y=y,
        hue=hue,
        palette=palette,
        # style=style,
        # markers=markers,
        edgecolor="none",
        legend=True,
        s=3,
        ax=ax
    )
    ax.axline((0, 0), slope=1, color="red", linewidth=1, linestyle="--")
    ax.set_aspect("equal")
    ax.set_xlim(0, 6000)
    ax.set_ylim(0, 6000)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks([0, 2000, 4000, 6000])
    ax.set_yticks([0, 2000, 4000, 6000])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_title(name)

show_legend = (0, fig_cols - 1)
for n_row in range(axes.shape[0]):
    for n_col in range(axes.shape[1]):
        ax = axes[n_row][n_col]
        if (n_row + 1) * (n_col + 1) > len(col_values):
            ax.remove()
            continue
        if (n_row, n_col) == show_legend:
            ax.legend(loc="center right", bbox_to_anchor=(1.5, 0.5))
        else:
            ax.legend_.remove()
        if n_row == (axes.shape[0] - 1):
            ax.set_xticklabels(["0", "2000", "4000", "6000"])
        if n_col == 0:
            ax.set_yticklabels(["0", "2000", "4000", "6000"])

fig.savefig("output/processed_read_length.png", dpi=900)


# In[ ]:




