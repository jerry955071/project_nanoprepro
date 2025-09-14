configfile: "config.json"

# include workflows
include: "workflows/Basecalling.smk"
include: "workflows/ReverseBasecalling.smk"
include: "workflows/NanoPreP.smk"
include: "workflows/Pychopper.smk"
include: "workflows/RATTLE.smk"
include: "workflows/BUSCO.smk"
include: "workflows/SQANTI.smk"
include: "workflows/Minimap2.smk"

# Custom functions
# 1. query: query from list of dict by key-value pair
from typing import List
def query(d:List[dict], k:str, v:str) -> dict:
    return [x for x in d if x[k] == v][0]

# 2. query_all: query from list of dict by key-value pair and return values from the requested field
def query_all(d:List[dict], k:str, v:str, k_out:str) -> List[str]:
    return [x[k_out] for x in d if x[k] == v]

# alias
docker_run="docker run --rm -u $(id -u) -v project_nanoprep_re:/mydata -w /mydata -v project_nanoprep:/mydata/old_data"

# wildcard constraints
wildcard_constraints:
    accuracy="sup|hac|fast|pre-called|bad|default|very-good",
    backend="phmm|edlib"



# ReverseBasecalling
rule ReverseBasecalling:
    input:
        expand(
            "outputs/ReverseBasecalling/longbow/{name}.json",
            name=[
                "mouse-testicle",
                "camellia-perpetua"
            ]
        )

# Run NanoPreP and Pychopper on simulated data 
# TODO: Rerun with NanoPreP-0.0.20
nanoprep_simulated = expand(
    "outputs/NanoPreP/mprof_nanoprep/beta{beta}/{name}_{accuracy}/full-length.fq",
    beta=["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0", "1.5", "2.0"],
    name=[
        "0.50-0.40-0.10",
        "0.50-0.45-0.05",
        "0.50-0.49-0.01",
        "0.50-0.50-0.00",
        "0.60-0.30-0.10",
        "0.60-0.35-0.05",
        "0.60-0.39-0.01",
        "0.60-0.40-0.00",
        "0.70-0.20-0.10",
        "0.70-0.25-0.05",
        "0.70-0.29-0.01",
        "0.70-0.30-0.00",
        "0.80-0.10-0.10",
        "0.80-0.15-0.05",
        "0.80-0.19-0.01",
        "0.80-0.20-0.00",
        "0.90-0.00-0.10",
        "0.90-0.05-0.05",
        "0.90-0.09-0.01",
        "0.90-0.10-0.00",
        "1.00-0.00-0.00"
    ],
    accuracy=["bad", "default", "very-good"]
)
pychopper_simulated = expand(
    "outputs/Pychopper/mprof_pychopper/{backend}/{name}_{accuracy}/full-length.fq",
    backend=["phmm", "edlib"],
    name=[
        "0.50-0.40-0.10",
        "0.50-0.45-0.05",
        "0.50-0.49-0.01",
        "0.50-0.50-0.00",
        "0.60-0.30-0.10",
        "0.60-0.35-0.05",
        "0.60-0.39-0.01",
        "0.60-0.40-0.00",
        "0.70-0.20-0.10",
        "0.70-0.25-0.05",
        "0.70-0.29-0.01",
        "0.70-0.30-0.00",
        "0.80-0.10-0.10",
        "0.80-0.15-0.05",
        "0.80-0.19-0.01",
        "0.80-0.20-0.00",
        "0.90-0.00-0.10",
        "0.90-0.05-0.05",
        "0.90-0.09-0.01",
        "0.90-0.10-0.00",
        "1.00-0.00-0.00"
    ],
    accuracy=["bad", "default", "very-good"]
)
rule simulated_data:
    input:
        nanoprep_simulated + pychopper_simulated
    shell:
        """
        echo "Simulated data done!"
        """

# Evaluation pipelines
evaluation_pipelines_input = expand(
    "outputs/{evaluation_pipeline}/{software}/{backend_or_beta}/{name}_{accuracy}",
    evaluation_pipeline=["BUSCO", "SQANTI"],
    software=["pychopper"],
    backend_or_beta=["phmm", "edlib"],
    name=[
        "ont-10x-human",
        "ont-visium-mouse",
        "ptr-109-bio1",
        "ptr-109-bio2",
        "ptr-111-bio1",
        "egr-109-bio1",
        "egr-109-bio2",
        "lch-109-bio1",
        "lch-109-bio2"
    ],
    accuracy=["sup", "hac", "fast"]
) + expand(
    "outputs/{evaluation_pipeline}/{software}/{backend_or_beta}/{name}_{accuracy}",
    evaluation_pipeline=["BUSCO", "SQANTI"],
    software=["nanoprep"],
    backend_or_beta=["beta0.1", "beta0.2", "beta0.3", "beta0.4", "beta0.5", "beta0.6", "beta0.7", "beta0.8", "beta0.9", "beta1.0", "beta1.5", "beta2.0"],
    name=[
        "ont-10x-human",
        "ont-visium-mouse",
        "ptr-109-bio1",
        "ptr-109-bio2",
        "ptr-111-bio1",
        "egr-109-bio1",
        "egr-109-bio2",
        "lch-109-bio1",
        "lch-109-bio2"
    ],
    accuracy=["sup", "hac", "fast"]
) + expand(
    "outputs/{evaluation_pipeline}/{software}/{backend_or_beta}/{name}_{accuracy}",
    evaluation_pipeline=["BUSCO", "SQANTI"],
    software=["nanoprep"],
    backend_or_beta=["beta0.1", "beta0.2", "beta0.3", "beta0.4", "beta0.5", "beta0.6", "beta0.7", "beta0.8", "beta0.9", "beta1.0", "beta1.5", "beta2.0"],
    name=[
        "mouse-retina-subset1",
        "mouse-retina-subset2"
    ],
    accuracy=["pre-called"]
)
evaluation_pipelines_input.remove("outputs/SQANTI/pychopper/edlib/ont-visium-mouse_fast") # failed to produce this output


# map processed/raw reads to genome
map2genome_input = expand(
    "outputs/Minimap2/minimap2_genome/{software}/{backend_or_beta}/{name}_{accuracy}.bam",
    software=["pychopper"],
    backend_or_beta=["phmm", "edlib"],
    name=[
        "ont-10x-human",
        "ont-visium-mouse",
        "ptr-109-bio1",
        "ptr-109-bio2",
        "ptr-111-bio1",
        "egr-109-bio1",
        "egr-109-bio2",
        "lch-109-bio1",
        "lch-109-bio2"
    ],
    accuracy=["sup", "hac", "fast"]
) + expand(
    "outputs/Minimap2/minimap2_genome/{software}/{backend_or_beta}/{name}_{accuracy}.bam",
    software=["nanoprep"],
    backend_or_beta=["beta0.1", "beta0.2", "beta0.3", "beta0.4", "beta0.5", "beta0.6", "beta0.7", "beta0.8", "beta0.9", "beta1.0", "beta1.5", "beta2.0"],
    name=[
        "ont-10x-human",
        "ont-visium-mouse",
        "ptr-109-bio1",
        "ptr-109-bio2",
        "ptr-111-bio1",
        "egr-109-bio1",
        "egr-109-bio2",
        "lch-109-bio1",
        "lch-109-bio2"
    ],
    accuracy=["sup", "hac", "fast"]
) + expand(
    "outputs/Minimap2/minimap2_genome/{software}/{backend_or_beta}/{name}_{accuracy}.bam",
    software=["nanoprep"],
    backend_or_beta=["beta0.1", "beta0.2", "beta0.3", "beta0.4", "beta0.5", "beta0.6", "beta0.7", "beta0.8", "beta0.9", "beta1.0", "beta1.5", "beta2.0"],
    name=[
        "mouse-retina-subset1",
        "mouse-retina-subset2"
    ],
    accuracy=["pre-called"]
) + expand(
    "outputs/Minimap2/minimap2_genome/raw/raw/{name}_{accuracy}.bam",
    name=[
        "mouse-retina-subset1",
        "mouse-retina-subset2"
    ],
    accuracy=["pre-called"]
) + expand(
    "outputs/Minimap2/minimap2_genome/raw/raw/{name}_{accuracy}.bam",
    name=[
        "ont-10x-human",
        "ont-visium-mouse",
        "ptr-109-bio1",
        "ptr-109-bio2",
        "ptr-111-bio1",
        "egr-109-bio1",
        "egr-109-bio2",
        "lch-109-bio1",
        "lch-109-bio2"
    ],
    accuracy=["sup", "hac", "fast"]
) 

rule evaluation_pipelines:
    input:
        evaluation_pipelines_input + map2genome_input

rule map2genome:
    input:
        map2genome_input

