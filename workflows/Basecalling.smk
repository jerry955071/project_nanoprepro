# Dorado 1.1.1
# Note: 
# 1. read splitting is allowed
def _model(wildcards):
    sample_dict = query(config["samples"], "name", wildcards.name)
    model = sample_dict["model"].replace("{accuracy}", wildcards.accuracy)
    return model

rule dorado:
    resources:
        gpu=1
    input:
        lambda wildcards: query(config["samples"], "name", wildcards.name)["signal_path"]
    output:
        directory("outputs/Basecalling/dorado/{name}/{accuracy}")
    params:
        model=_model
    log:
        "logs/Basecalling/dorado/{name}/{accuracy}.log"
    shell:
        """
        docker --context docker-olaf \
            run --runtime=nvidia --gpus all -u $(id -u) --rm -v project_nanoprep_re:/mydata -w /mydata nanoporetech/dorado:shae423e761540b9d08b526a1eb32faf498f32e8f22 \
                dorado basecaller \
                    -x cuda:all \
                    -r \
                    -v \
                    --emit-fastq \
                    -o {output} \
                    --no-trim \
                    {params.model} {input} \
                    2> {log} \
                    1> {log}
        """

# Dorado 0.9.5: nanoporetech/dorado:sha268dcb4cd02093e75cdc58821f8b93719c4255ed
rule dorado_legacy:
    resources:
        gpu=1
    input:
        lambda wildcards: query(config["samples"], "name", wildcards.name)["signal_path"]
    output:
        directory("outputs/Basecalling/dorado-legacy/{name}/{accuracy}")
    params:
        model=lambda wildcards: {
            "sup": "/models/dna_r9.4.1_e8_sup@v3.6",
            "hac": "/models/dna_r9.4.1_e8_hac@v3.3",
            "fast": "/models/dna_r9.4.1_e8_fast@v3.4"
        }[wildcards.accuracy]
    log:
        "logs/Basecalling/dorado/{name}/{accuracy}.log"
    shell:
        """
        docker --context docker-olaf \
            run --runtime=nvidia --gpus all -u $(id -u) --rm -v project_nanoprep_re:/mydata -w /mydata \
            nanoporetech/dorado:sha268dcb4cd02093e75cdc58821f8b93719c4255ed \
                dorado basecaller \
                    -x cuda:all \
                    -r \
                    -v \
                    --emit-fastq \
                    -o {output} \
                    --no-trim \
                    {params.model} {input} \
                    2> {log} \
                    1> {log}
        """

# merge basecalling results to single FASTQ file
rule aggregate_fastq:
    input:
        lambda wildcards: {
            "R10.4.1": f"outputs/Basecalling/dorado/{wildcards.name}/{wildcards.accuracy}",
            "R9.4.1": f"outputs/Basecalling/dorado-legacy/{wildcards.name}/{wildcards.accuracy}"
        }[query(config["samples"], "name", wildcards.name)["chemistry"]]
    output:
        "outputs/Basecalling/aggregate_fastq/{name}_{accuracy}.fq"
    log:
        "logs/Basecalling/aggregate_fastq/{name}_{accuracy}.log"
    shell:
        """
        echo "Listing all FASTQ files:" > {log}
        find {input} -name "*.fastq" >> {log}
        find {input} -name "*.fastq" -exec cat {{}} \; > {output}
        """

