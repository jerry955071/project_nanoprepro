def _full_length(wildcards):
    if wildcards.software == "nanoprep":
        return "outputs/NanoPreP/mprof_nanoprep/{backend_or_beta}/{name}_{accuracy}/full-length.fq"
    elif wildcards.software == "pychopper":
        return "outputs/Pychopper/mprof_pychopper/{backend_or_beta}/{name}_{accuracy}/full-length.fq"
    else:
        raise ValueError("Unsupported software: {}".format(wildcards.software))

rule cluster:
    threads: 8
    input:
        _full_length
    output:
        fout="outputs/RATTLE/cluster/{software}/{backend_or_beta}/{name}_{accuracy}/clusters.out",
        outdir=directory("outputs/RATTLE/cluster/{software}/{backend_or_beta}/{name}_{accuracy}/")
    params:
        "--iso --rna"
    log:
        "logs/RATTLE/cluster/{software}/{backend_or_beta}/{name}_{accuracy}.log"
    shell:
        """
        mkdir -p {output.outdir}
        {docker_run} chiaenu/rattle:1.0 \
            cluster \
                {params} \
                -i {input} \
                -t {threads} \
                -o {output.outdir} \
            2> {log} \
            1> {log} 
        """

rule correct:
    threads: 8
    input:
        cluster="outputs/RATTLE/cluster/{software}/{backend_or_beta}/{name}_{accuracy}/clusters.out",
        filtered=_full_length
    output:
        outdir=directory("outputs/RATTLE/correct/{software}/{backend_or_beta}/{name}_{accuracy}"),
        corr="outputs/RATTLE/correct/{software}/{backend_or_beta}/{name}_{accuracy}/corrected.fq",
        ucrr="outputs/RATTLE/correct/{software}/{backend_or_beta}/{name}_{accuracy}/uncorrected.fq",
        cons="outputs/RATTLE/correct/{software}/{backend_or_beta}/{name}_{accuracy}/consensi.fq"
    log:
        "logs/RATTLE/correct/{software}/{backend_or_beta}/{name}_{accuracy}.log"
    shell:
        """
        mkdir -p {output.outdir}
        {docker_run} chiaenu/rattle:1.0 \
             correct \
                -i {input.filtered} \
                -c {input.cluster} \
                -o {output.outdir} \
                -t {threads} \
            2> {log} \
            1> {log} 
        """

rule polish:
    threads: 8
    input:
        cons="outputs/RATTLE/correct/{software}/{backend_or_beta}/{name}_{accuracy}/consensi.fq"
    output:
        fout="outputs/RATTLE/polish/{software}/{backend_or_beta}/{name}_{accuracy}/transcriptome.fq",
        outdir=directory("outputs/RATTLE/polish/{software}/{backend_or_beta}/{name}_{accuracy}")
    log:
        "logs/RATTLE/polish/{software}/{backend_or_beta}/{name}_{accuracy}.log"
    shell:
        """
        mkdir -p {output.outdir}
        {docker_run} chiaenu/rattle:1.0 \
            polish \
                -i {input.cons} \
                -o {output.outdir} \
                --rna \
                -t {threads} \
            2> {log} \
            1> {log} 
        """