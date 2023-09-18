# rules to initialize database
localrules: init, emu_prebuilt, download_silva, filt_silva, silva2ncbi_map

def init_db(db=CLASSIFIER):
    if db == "emu":
        return expand(DATABASE_DIR + "/emu/" + DATABASE_PREBUILT + "_prebuilt/{file}", file = ["species_taxid.fasta", "taxonomy.tsv"])
    elif db == "minimap2lca":
        return expand(DATABASE_DIR + "/silva/{file}", file = ["SILVA_ssu_nr99_filt.mmi", "silva_to_ncbi.map"])
    elif db == "blast2lca":
        return multiext(DATABASE_DIR + "/silva/SILVA_ssu_nr99_filt.fasta", 
            ".ndb",
            ".nhr",
            ".nin",
            ".not",
            ".nsq",
            ".ntf",
            ".nto"
            ), DATABASE_DIR + "/silva/silva_to_ncbi.map"
    else:
        raise ValueError("Classifier not supported")

rule init:
    input: init_db() 

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

def add_spikein(fasta_in, tax_in, spikein, fasta_out, tax_out, spikein_tax):
    # open spikein fasta and rename the header
    # add spikein to the end of the prebuilt fasta, and make new unique headers 
    # open spikein taxonomy and add to the end of the prebuilt taxonomy
    # ugly but works: 1111111111
    # custom spikein_tax according to estimation appraoch
    with open(fasta_out, "w") as f1, open(tax_out, "w") as f2:
        with open(fasta_in, "r") as g1:
            for line in g1:
                f1.write(line)
        with open(tax_in, "r") as g2:
            for line in g2:
                f2.write(line)
        with open(spikein, "r") as g:
            i = 0
            suffix = 1111111111
            for line in g:
                if line.startswith(">"):
                    line = ">" + str(suffix + i) + "\n"
                    line2 = str(suffix + i) + "\t" + str(spikein_tax) + "\n" 
                    f2.write(line2)
                    i += 1 
                f1.write(line)

rule emu_spikein:
    input: 
        rules.emu_prebuilt.output,
    output:
        expand(DATABASE_DIR + "/emu/" + DATABASE_PREBUILT + "_prebuilt_spikein/{file}", file = ["species_taxid.fasta", "taxonomy.tsv"])
    params:
        taxdb = DATABASE_PREBUILT,
        database_dir = DATABASE_DIR,
        spikein = get_spikein_fasta(),
    message: "Adding spikein to the prebuit Emu database [{params.taxdb}]"
    run:
        add_spikein(input[0], input[1], params.spikein, output[0], output[1], spikein_tax="spikein")
       
# download silva database
rule download_silva:
    output:
        temp(expand(DATABASE_DIR + "/silva/{file}", file = ["SILVA_ssu_nr99_full.fasta", "tax_ncbi-species_nr99.txt", "taxmap_ssu_ref_nr99.txt"]))
    message: "Downloading the Silva database"
    conda: "../envs/lca.yaml"
    params:
        database_dir = DATABASE_DIR,
        silvaFastaURL = "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz",
        silvaTaxNcbiSpURL = "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/ncbi/tax_ncbi-species_ssu_ref_nr99_138.1.txt.gz",
        silvaTaxmapURL = "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_ssu_ref_nr_138.1.txt.gz",
    log: "logs/taxonomy/silva/database.log"
    benchmark: "benchmarks/taxonomy/silva/database.txt"
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["simple"],
    shell: 
        """
        mkdir -p {params.database_dir}/silva
        wget -O {params.database_dir}/silva/SILVA_ssu_nr99_full.fasta.gz {params.silvaFastaURL} 1> {log} 2>&1
        wget -O {params.database_dir}/silva/tax_ncbi-species_nr99.txt.gz {params.silvaTaxNcbiSpURL} 1>> {log} 2>&1
        wget -O {params.database_dir}/silva/taxmap_ssu_ref_nr99.txt.gz {params.silvaTaxmapURL} 1>> {log} 2>&1
        gunzip {params.database_dir}/silva/SILVA_ssu_nr99_full.fasta.gz 1>> {log} 2>&1
        gunzip {params.database_dir}/silva/tax_ncbi-species_nr99.txt.gz 1>> {log} 2>&1
        gunzip {params.database_dir}/silva/taxmap_ssu_ref_nr99.txt.gz 1>> {log} 2>&1
        """

# filt silva database;
# no Phage & Eukaryota & Mitochondria & Chloroplast
# no Calanus (wrong annotation in SILVA)
# rna -> dna

rule filt_silva:
    input: DATABASE_DIR + "/silva/SILVA_ssu_nr99_full.fasta"
    output: DATABASE_DIR + "/silva/SILVA_ssu_nr99_filt.fasta"
    message: "Filtering the Silva database"
    conda: "../envs/lca.yaml"
    log: "logs/taxonomy/silva/filt.log"
    benchmark: "benchmarks/taxonomy/silva/filt.txt"
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["simple"],
    shell: 
        """
        seqkit grep -v -r -i -n -p Phage -p Eukaryota -p Mitochondria -p Chloroplast -p Calanus {input} | seqkit seq --rna2dna -o {output} 1> {log} 2>&1
        """

rule silva2ncbi_map:
    input: 
        ncbisp = DATABASE_DIR + "/silva/tax_ncbi-species_nr99.txt",
        taxmap = DATABASE_DIR + "/silva/taxmap_ssu_ref_nr99.txt"
    output: DATABASE_DIR + "/silva/silva_to_ncbi.map"
    message: "Generating the Silva to NCBI synonyms map"
    conda: "../envs/lca.yaml"
    log: "logs/taxonomy/silva/silva_to_ncbi_map.log"
    benchmark: "benchmarks/taxonomy/silva/silva_to_ncbi_map.txt"
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["simple"],
    shell: 
        """
        Rscript {workflow.basedir}/scripts/generateSynonyms.R \
		--ncbisp {input.ncbisp} \
		--taxmap  {input.taxmap} \
		--out {output} 1> {log} 2>&1
        """

rule silva_spikein:
    input:
        rules.filt_silva.output,
        rules.silva2ncbi_map.output,
    output:
        expand(DATABASE_DIR + "/silva_spikein/{file}", file = ["SILVA_ssu_nr99_filt.fasta", "silva_to_ncbi.map"])
    params:
        spikein = get_spikein_fasta(),
    message: "Adding spikein to the Silva database"
    run:
        # similar to rules.emu_spikein
        # https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=10710&lvl=3&lin=f&keep=1&srchmode=1&unlock
        # use lambda phage taxid:10710
        add_spikein(input[0], input[1], input.spikein, output[0], output[1], spikein_tax="10710")

rule minimap2silva_db:
    input: DATABASE_DIR + "/{silva}/SILVA_ssu_nr99_filt.fasta"
    output: DATABASE_DIR + "/{silva}/SILVA_ssu_nr99_filt.mmi"
    message: "Building the minimap2 database with SILVA reference"
    conda: "../envs/lca.yaml"
    params:
        k = 15,
    log: "logs/taxonomy/silva/minimap2{silva}_db.log"
    benchmark: "benchmarks/taxonomy/silva/minimap2{silva}_db.txt"
    threads: config["threads"]["large"]
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["simple"],
    shell: 
        """
        minimap2 -d {output} -t {threads} -k {params.k} {input} 1> {log} 2>&1 
        """

rule blast2silva_db:
    input: DATABASE_DIR + "/{silva}/SILVA_ssu_nr99_filt.fasta"
    output: 
        multiext(DATABASE_DIR + "/{silva}/SILVA_ssu_nr99_filt.fasta",
            ".ndb",
            ".nhr",
            ".nin",
            ".not",
            ".nsq",
            ".ntf",
            ".nto"
            )
    message: "Building the blast database with SILVA reference"
    conda: "../envs/lca.yaml"
    log: "logs/taxonomy/silva/blast2{silva}_db.log"
    benchmark: "benchmarks/taxonomy/silva/blast2{silva}_db.txt"
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["simple"],
    shell: "makeblastdb -in {input} -dbtype nucl 1> {log} 2>&1" 
