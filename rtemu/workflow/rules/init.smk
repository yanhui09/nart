# rules to initialize database
localrules: init_database, emu_prebuilt
rule init:
    input: 
        expand(DATABASE_DIR + "/emu/" + DATABASE_PREBUILT + "_prebuilt/{file}", file = ["species_taxid.fasta", "taxonomy.tsv"])

# use prebuilt Emu database
rule emu_prebuilt:
    output:
        expand(DATABASE_DIR + "/emu/" + DATABASE_PREBUILT + "_prebuilt/{file}", file = ["species_taxid.fasta", "taxonomy.tsv"])
    message: "Downloading the prebuit Emu database [{params.taxdb}]"
    conda: "../envs/emu.yaml"
    params:
        taxdb = DATABASE_PREBUILT,
        database_dir = DATABASE_DIR,
    log: "logs/taxonomy/emu/database_{paramas.taxdb}.log"
    benchmark: "benchmarks/taxonomy/emu/databases_{params.taxdb}.txt"
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["simple"],
    shell: 
        """
        mkdir -p {params.database_dir}/emu/{params.taxdb}_prebuilt
        osf -p 56uf7 fetch osfstorage/emu-prebuilt/{params.taxdb}.tar {params.database_dir}/emu/{params.taxdb}_prebuilt/{params.taxdb}.tar 1> {log} 2>&1
        tar -xvf {params.database_dir}/emu/{params.taxdb}_prebuilt/{params.taxdb}.tar -C {params.database_dir}/emu/{params.taxdb}_prebuilt 1>> {log} 2>&1
        rm {params.database_dir}/emu/{params.taxdb}_prebuilt/{params.taxdb}.tar
        """
