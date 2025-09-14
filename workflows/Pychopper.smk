def _fastq(wildcards):
    # get sample dict
    sample_dict = query(d=config["samples"], k="name", v=wildcards.name)
    
    # get fastq path if species == "ptr-simulated"
    if sample_dict.get("species") == "ptr-simulated":
        return sample_dict["fastq_path"].replace("{accuracy}", wildcards.accuracy)

    # get fastq path from config if exists
    if sample_dict.get("fastq_path"):
        return sample_dict["fastq_path"]

    # otherwise, use basecalled fastq create by workflow/Basecalling.smk
    if sample_dict.get("signal_path"):
        return f"outputs/Basecalling/aggregate_fastq/{wildcards.name}_{wildcards.accuracy}.fq"

    return None

rule mprof_pychopper:
    threads: 1
    input:
        _fastq
    output:
        fl="outputs/Pychopper/mprof_pychopper/{backend}/{name}_{accuracy}/full-length.fq",
        fu="outputs/Pychopper/mprof_pychopper/{backend}/{name}_{accuracy}/chimeric.fq",
        tr="outputs/Pychopper/mprof_pychopper/{backend}/{name}_{accuracy}/truncated.fq",
        report="outputs/Pychopper/mprof_pychopper/{backend}/{name}_{accuracy}/report.pdf",
        stats="outputs/Pychopper/mprof_pychopper/{backend}/{name}_{accuracy}/stats.tsv",
        mprof="outputs/Pychopper/mprof_pychopper/{backend}/{name}_{accuracy}/mprof.txt"
    params:
        kit=lambda wildcards: query(d=config["samples"], k="name", v=wildcards.name)["kit"],
        min_qual=7,
        min_len=0,
        backend=lambda wildcards: wildcards.backend,
        autotune_nr=100000,
        rescue=lambda wildcards: query(d=config["samples"], k="name", v=wildcards.name)["kit"]
    log:
        "logs/Pychopper/mprof_nanoprep/{backend}/{name}_{accuracy}.log"
    shell:
        """
        {docker_run} chiaenu/pychopper:2.7.10 \
            mprof run \
                --include-children \
                --multiprocess \
                --output {output.mprof} \
                pychopper \
                    -k {params.kit} \
                    -Q {params.min_qual} \
                    -z {params.min_len} \
                    -r {output.report} \
                    -u {output.tr} \
                    -w {output.fu} \
                    -S {output.stats} \
                    -Y {params.autotune_nr} \
                    -m {params.backend} \
                    -x {params.rescue} \
                    -t {threads} \
                    {input} \
                    {output.fl} \
            2> {log} \
            1> {log}
        """
    

rule cutadapt:
    threads: 4
    input:
        "outputs/Pychopper/mprof_pychopper/{backend}/{name}_{accuracy}/full-length.fq"
    output:
        "outputs/Pychopper/cutadapt/{backend}/{name}_{accuracy}/full-length.fq"
    log:
        "logs/Pychopper/cutadapt/{backend}/{name}_{accuracy}.log"
    shell:
        """
        {docker_run} chiaenu/cutadapt:5.1 \
            cutadapt \
                --poly-a \
                -o {output} \
                {input} \
            2> {log} \
            1> {log} 
        """