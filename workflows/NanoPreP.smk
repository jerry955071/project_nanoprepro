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

rule mprof_nanoprep:
    threads: 1
    input:
        _fastq
    output:
        fl="outputs/NanoPreP/mprof_nanoprep/beta{beta}/{name}_{accuracy}/full-length.fq",
        ch="outputs/NanoPreP/mprof_nanoprep/beta{beta}/{name}_{accuracy}/chimeric.fq",
        tr="outputs/NanoPreP/mprof_nanoprep/beta{beta}/{name}_{accuracy}/truncated.fq",
        report="outputs/NanoPreP/mprof_nanoprep/beta{beta}/{name}_{accuracy}/report.html",
        mprof="outputs/NanoPreP/mprof_nanoprep/beta{beta}/{name}_{accuracy}/mprof.txt"
    params:
        p5_sense=lambda wildcards: query(d=config["samples"], k="name", v=wildcards.name)["p5_sense"],
        p3_sense=lambda wildcards: query(d=config["samples"], k="name", v=wildcards.name)["p3_sense"],
        skip_lowq=0,
        filter_lowq=7,
        filter_short=0,
        beta=lambda wildcards: wildcards.beta,
        n=100000,
        flags="--trim_adapter --trim_poly --poly_w 6 --poly_k 4 --pid_body 0.8",
        orientation=1
    log:
        "logs/NanoPreP/mprof_nanoprep/beta{beta}/{name}_{accuracy}.log"
    shell:
        """
        {docker_run} nanoprep-mprof:0.0.20 \
            mprof run \
                --include-children \
                --multiprocess \
                --output {output.mprof} \
                nanoprep \
                    --beta {params.beta} \
                    --p5_sense {params.p5_sense} \
                    --p3_sense {params.p3_sense} \
                    --skip_lowq {params.skip_lowq} \
                    --filter_lowq {params.filter_lowq} \
                    --filter_short {params.filter_short} \
                    --input_fq {input} \
                    --output_full_length {output.fl} \
                    --output_fusion {output.ch} \
                    --output_truncated {output.tr} \
                    --report {output.report} \
                    --processes {threads} \
                    -n {params.n} \
                    {params.flags} \
            2> {log} \
            1> {log}
        """
