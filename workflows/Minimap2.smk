def _full_length(wildcards):
    if wildcards.software == "nanoprep":
        return f"outputs/NanoPreP/mprof_nanoprep/{wildcards.backend_or_beta}/{wildcards.name}_{wildcards.accuracy}/full-length.fq"
    elif wildcards.software == "pychopper":
        return f"outputs/Pychopper/mprof_pychopper/{wildcards.backend_or_beta}/{wildcards.name}_{wildcards.accuracy}/full-length.fq"
    elif wildcards.software == "raw":
        if wildcards.accuracy == "pre-called":
            return query(d=config["samples"], k="name", v=wildcards.name)["fastq_path"]
        else:
            return f"outputs/Basecalling/aggregate_fastq/{wildcards.name}_{wildcards.accuracy}.fq"
    else:
        raise ValueError("Unsupported software: {}".format(wildcards.software))

# path to minimap2 executable
minimap2 = "src/minimap2-2.30/minimap2-2.30_x64-linux/minimap2"

# map processed reads to genome
rule minimap2_genome:
    threads: 8
    input:
        _full_length
    output:
        "outputs/Minimap2/minimap2_genome/{software}/{backend_or_beta}/{name}_{accuracy}.bam"
    params:
        genome=lambda wildcards: config["references"][query(config["samples"], "name", wildcards.name)["species"]]["fasta"]
    log:
        "logs/Minimap2/minimap2_genome/{software}/{backend_or_beta}/{name}_{accuracy}.log"
    shell:
        """
        {minimap2} \
            -ax splice --secondary=no --eqx \
            -t {threads} \
            -o {output} \
            {params.genome} \
            {input} \
            2> {log} 1> {log}
        """
