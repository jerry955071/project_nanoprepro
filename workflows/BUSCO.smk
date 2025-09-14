def _lineage_dataset(wildcards):
    sample_dict = query(config["samples"], "name", wildcards.name)
    lineages = {
        "egr": "eudicotyledons_odb12",
        "ptr": "malpighiales_odb12",
        "lch": "embryophyta_odb12",
        "hsa": "primates_odb12",
        "mmu": "rodentia_odb12",
    }
    if lineages.get(sample_dict["species"]):
        return f"outputs/BUSCO/lineage_data/{lineages[sample_dict['species']]}"
    else:
        raise ValueError(f"Lineage dataset for species {sample_dict['species']} not found.")

rule busco:
    threads: 10
    input:
        fq="outputs/RATTLE/polish/{software}/{backend_or_beta}/{name}_{accuracy}/transcriptome.fq",
        lineage_dataset=_lineage_dataset
    output:
        fa=temp("outputs/BUSCO/{software}/{backend_or_beta}/{name}_{accuracy}.fa"),
        out_dir=directory("outputs/BUSCO/{software}/{backend_or_beta}/{name}_{accuracy}/")
    log:
        "logs/BUSCO/busco/{software}/{backend_or_beta}/{name}_{accuracy}.log"
    shell:
        """
        # Convert fq to fa
        sed -n '1~4s/^@/>/p;2~4p' {input.fq} > {output.fa}
        
        # Run BUSCO
        {docker_run} ezlabgva/busco:v6.0.0_cv1 \
            busco \
                -m transcriptome \
                -i {output.fa} \
                -l {input.lineage_dataset} \
                -o {output.out_dir} \
                --cpu {threads} \
                --offline \
            2> {log} \
            1> {log}
        """

rule download_lineage_data:
    threads: 8
    output:
        directory("outputs/BUSCO/lineage_data/{lineage}")
    params:
        lineage=lambda wildcards: {
            "eudicotyledons_odb12": "https://busco-data.ezlab.org/v5/data/lineages/eudicotyledons_odb12.2025-07-01.tar.gz",
            "malpighiales_odb12": "https://busco-data.ezlab.org/v5/data/lineages/malpighiales_odb12.2025-07-01.tar.gz",
            "embryophyta_odb12": "https://busco-data.ezlab.org/v5/data/lineages/embryophyta_odb12.2025-07-01.tar.gz",
            "primates_odb12": "https://busco-data.ezlab.org/v5/data/lineages/primates_odb12.2025-07-01.tar.gz",
            "rodentia_odb12": "https://busco-data.ezlab.org/v5/data/lineages/rodentia_odb12.2025-07-01.tar.gz",
        }[wildcards.lineage]
    log:
        "logs/BUSCO/download_lineage_data/{lineage}.log"
    shell:
        """
        wget -O - -o {log} {params.lineage} | tar -xz -C outputs/BUSCO/lineage_data/ >> {log} 2>&1
        """