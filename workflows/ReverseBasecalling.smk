rule longbow:
    threads: 24
    input:
        lambda wildcards: query(config["samples"], "name", wildcards.name)["fastq_path"]
    output:
        "outputs/ReverseBasecalling/longbow/{name}.json"
    log: 
        "logs/ReverseBasecalling/longbow/{name}.log"
    shell:
        """
        docker run --rm -v project_nanoprep_re:/mydata -w /mydata chiaenu/longbow:2.3.1 \
            longbow -t {threads} -i {input} -o {output} -b \
            2> {log} \
            1> {log}

        docker run --rm -v project_nanoprep_re:/mydata -w /mydata chiaenu/longbow:2.3.1 \
            chown $(id -u):$(id -g) {output}
        """