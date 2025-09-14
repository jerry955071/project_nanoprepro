rule sqanti:
    threads: 8
    input:
        isoforms="outputs/RATTLE/polish/{software}/{backend_or_beta}/{name}_{accuracy}/transcriptome.fq",
        refGTF=lambda wildcards: config["references"][query(config["samples"], "name", wildcards.name)["species"]]["gtf"],
        refFasta=lambda wildcards: config["references"][query(config["samples"], "name", wildcards.name)["species"]]["fasta"]
    output:
        out_dir=directory("outputs/SQANTI/{software}/{backend_or_beta}/{name}_{accuracy}/")
    log:
        "logs/SQANTI/sqanti/{software}/{backend_or_beta}/{name}_{accuracy}.log"
    shell:
        """
        {docker_run} anaconesalab/sqanti3:v5.5.1 \
            sqanti3_qc.py \
                --isoforms {input.isoforms} \
                --refGTF {input.refGTF} \
                --refFasta {input.refFasta} \
                --fasta \
                --force_id_ignore \
                --aligner_choice minimap2 \
                -d {output.out_dir} \
                -t {threads} \
            2> {log} \
            1> {log}
        """
        