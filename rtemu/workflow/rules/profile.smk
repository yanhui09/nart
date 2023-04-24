rule emu:
    input: 
        rules.emu_prebuilt.output,
        fq = rules.q_filter.output,
    output: temp("{batch}/emu/{barcode}_rel-abundance.tsv")
    conda: "../envs/emu.yaml"
    params: 
        db = DATABASE_DIR + "/emu/" + DATABASE_PREBUILT + "_prebuilt",
        outdir = os.getcwd() + "/{batch}/emu"
    log: "logs/emu/{barcode}_{batch}.log"
    benchmark: "benchmarks/emu/{barcode}_{batch}.txt"
    threads: config["threads"]["large"]
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["default"],
    shell:
        "emu abundance --db {params.db} {input.fq} "
        " --keep-counts --threads {threads} "
        " --output-dir {params.outdir} "
        " > {log} 2>&1 "

def get_profiles(wildcards, classifier=CLASSIFIER):
    barcodes = get_qced_barcodes(wildcards) 
    if classifier == "emu":
        return expand("{{batch}}/emu/{barcode}_rel-abundance.tsv", barcode=barcodes)
    else:
        return expand("{{batch}}/{{classifier}}/{barcode}.ri", barcode=barcodes)

localrules: emu_merge, lca_merge
rule emu_merge:
    input: 
        demux_dir = "{batch}/demultiplexed",
        otutab = lambda wildcards: get_profiles(wildcards, classifier="emu"),
    output: "batches/emu_{batch}.tsv"
    params:
        db = DATABASE_PREBUILT,
    resources:
        mem = config["mem"]["normal"],
    run:
        import pandas as pd
        # merge tsv files
        otu_table = pd.DataFrame()
        for f in input.otutab:
            # get file name
            barcode = os.path.basename(f).split("_")[0]
            table = pd.read_csv(f, sep="\t")
            if params.db == 'silva':
                # siliva ouput one 'lineage' column
                table = table.rename(columns={"lineage": "taxonomy"})
            else:
                # combine "superkingdom", "phylum", "class", "order", "family", "genus", "species" into "taxonomy"
                # if empty, skip
                columns_to_combine = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
                # na to ""
                table[columns_to_combine] = table[columns_to_combine].fillna("")
                table["taxonomy"] = table[columns_to_combine].apply(lambda x: ";".join(x), axis=1)
            # if tax_id is "unassigned", taxonomy is "unassigned"
            table.loc[table["tax_id"] == "unassigned", "taxonomy"] = "unassigned"
            # only keep "tax_id", "estimated counts", "taxonomy"
            table = table[["tax_id", "taxonomy", "estimated counts"]]
            # use integer for "estimated counts"
            table["estimated counts"] = table["estimated counts"].astype(int)
            # rename "estimated counts" to sample name
            table = table.rename(columns={"estimated counts": barcode})
            # merge table into otu_table
            if otu_table.empty:
                otu_table = table
            else:
                otu_table = pd.merge(otu_table, table, on=["tax_id", "taxonomy"], how="outer")
        # fill NaN with 0
        otu_table = otu_table.fillna(0)
        # write otu_table to file; integer without decimal
        otu_table.to_csv(output[0], sep="\t", index=False, float_format="%.0f")

# minimap2lca
rule minimap2:
    input:
        mmi = rules.minimap2silva_db.output,
        fq = rules.q_filter.output,
    output: temp("{batch}/minimap2lca/{barcode}.sam")
    conda: "../envs/lca.yaml"
    params:
        K = '500M',
        x = 'map-ont',
        f = 10000,
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["default"],
    log: "logs/minimap2lca/{barcode}_{batch}_minimap2.log"
    benchmark: "benchmarks/minimap2lca/{barcode}_{batch}_minimap2.txt"
    threads: config["threads"]["large"]
    shell:
        """
        minimap2 -t {threads} -K {params.K} -ax {params.x} -f {params.f} --secondary=no {input.mmi} {input.fq} > {output} 2> {log}
        """

rule minimap2rma:
    input:
        fq = rules.q_filter.output,
        sam = rules.minimap2.output,
        silva2ncbi_map = rules.silva2ncbi_map.output,
    output: temp("{batch}/minimap2lca/{barcode}.rma")
    conda: "../envs/lca.yaml"
    params:
        maxMatchesPerRead = config["lca"]["maxMatchesPerRead"],
        topPercent = config["lca"]["topPercent"],
        minSupportPercent = config["lca"]["minSupportPercent"],
        minPercentReadCover = config["lca"]["minPercentReadCover"],
        minPercentReferenceCover = config["lca"]["minPercentReferenceCover"],
        lcaAlgorithm = config["lca"]["lcaAlgorithm"],
        lcaCoveragePercent = config["lca"]["lcaCoveragePercent"],
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["default"],
    log: "logs/minimap2lca/{barcode}_{batch}_sam2rma.log"
    benchmark: "benchmarks/minimap2lca/{barcode}_{batch}_sam2rma.txt"
    threads: config["threads"]["large"]
    shell:
        """
        sam2rma \
        -i {input.sam} \
        -r {input.fq} \
        -o {output} \
        -lg \
        -m {params.maxMatchesPerRead} \
        -top {params.topPercent} \
        -supp {params.minSupportPercent} \
        -mrc {params.minPercentReadCover} \
        -mrefc {params.minPercentReferenceCover} \
        -alg {params.lcaAlgorithm} \
        -lcp {params.lcaCoveragePercent} \
        -ram readCount \
        -s2t {input.silva2ncbi_map} > {log} 2>&1
        """

rule blastn:
    input:
        multiext(DATABASE_DIR + "/silva/SILVA_ssu_nr99_filt.fasta",
            ".ndb",
            ".nhr",
            ".nin",
            ".not",
            ".nsq",
            ".ntf",
            ".nto"
            ),
        db = rules.filt_silva.output,
        fq = rules.q_filter.output,
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["default"],
    output:
        fa = temp("{batch}/blast2lca/{barcode}.fasta"), 
        tab = temp("{batch}/blast2lca/{barcode}.tab")
    conda: "../envs/lca.yaml"
    params:
        max_target_seqs = 10,
    log: "logs/blast2lca/{barcode}_{batch}.log"
    benchmark: "benchmarks/blast2lca/{barcode}_{batch}.txt"
    threads: config["threads"]["large"]
    shell: 
        """
        seqkit fq2fa {input.fq} -o {output.fa} > {log} 2>&1
        blastn -query {output.fa} -db {input.db} -out {output.tab} -outfmt 6 -num_threads {threads} -max_target_seqs {params.max_target_seqs} >> {log} 2>&1
        """

rule blast2rma:
    input:
        fa = rules.blastn.output.fa,
        tab = rules.blastn.output.tab,
        silva2ncbi_map = rules.silva2ncbi_map.output,
    output: temp("{batch}/blast2lca/{barcode}.rma")
    conda: "../envs/lca.yaml"
    params:
        maxMatchesPerRead = config["lca"]["maxMatchesPerRead"],
        topPercent = config["lca"]["topPercent"],
        minSupportPercent = config["lca"]["minSupportPercent"],
        minPercentReadCover = config["lca"]["minPercentReadCover"],
        minPercentReferenceCover = config["lca"]["minPercentReferenceCover"],
        lcaAlgorithm = config["lca"]["lcaAlgorithm"],
        lcaCoveragePercent = config["lca"]["lcaCoveragePercent"],
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["default"],
    log: "logs/blast2lca/{barcode}_{batch}_blast2rma.log"
    benchmark: "benchmarks/blast2lca/{barcode}_{batch}_blast2rma.txt"
    threads: config["threads"]["large"]
    shell:
        """
        blast2rma \
        -i {input.tab} \
        -f BlastTab \
        -bm BlastN \
        -r {input.fa} \
        -o {output} \
        -lg \
        -m {params.maxMatchesPerRead} \
        -top {params.topPercent} \
        -supp {params.minSupportPercent} \
        -mrc {params.minPercentReadCover} \
        -mrefc {params.minPercentReferenceCover} \
        -alg {params.lcaAlgorithm} \
        -lcp {params.lcaCoveragePercent} \
        -ram readCount \
        -s2t {input.silva2ncbi_map} > {log} 2>&1
        """

#rule get_read_info:
#    input: rules.megan_lca.output
#    output: 
#        ncbi = temp("{batch}/minimap2lca/{barcode}.ncbi"),
#        path = temp("{batch}/minimap2lca/{barcode}.path"),
#        ri = temp("{batch}/minimap2lca/{barcode}.ri"),
#    conda: "../envs/lca.yaml"
#    resources:
#        mem = config["mem"]["normal"],
#        time = config["runtime"]["default"],
#    log: "logs/minimap2lca/{barcode}_{batch}_readinfo.log"
#    benchmark: "benchmarks/minimap2lca/{barcode}_{batch}_readinfo.txt"
#    shell:
#        """
#        rma2info \
#        -i {input} \
#        -r2c Taxonomy \
#        -n -r -mro \
#        -o {output.path} > {log} 2>&1
#        rma2info \
#        -i {input} \
#        -r2c Taxonomy \
#        -p -mro \
#        -o {output.ncbi} >> {log} 2>&1
#        join -t $'\\t' \
#        <(sort {output.ncbi}) \
#        <(sort {output.path}) \
#        > {output.ri}
#        """
wildcard_constraints:
    classifier = '|'.join(["minimap2lca", "blast2lca"])

rule rma2counts:
    input: "{batch}/{classifier}/{barcode}.rma"
    output: 
        ncbi = temp("{batch}/{classifier}/{barcode}.ncbi"),
        path = temp("{batch}/{classifier}/{barcode}.path"),
        ri = temp("{batch}/{classifier}/{barcode}.ri"),
    conda: "../envs/lca.yaml"
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["default"],
    log: "logs/{classifier}/{barcode}_{batch}_readinfo.log"
    benchmark: "benchmarks/{classifier}/{barcode}_{batch}_readinfo.txt"
    shell:
        """
        rma2info \
        -i {input} \
        -c2c Taxonomy \
        -mro \
        -o {output.path} > {log} 2>&1
        rma2info \
        -i {input} \
        -c2c Taxonomy \
        -p -mro \
        -o {output.ncbi} >> {log} 2>&1
        paste \
        <(cut -f1 {output.path}) \
        <(sed 's/\\[[A-Z]\\] //g;s/; /;/g' {output.ncbi}) \
        > {output.ri}
        """

rule lca_merge:
    input: 
        demux_dir = "{batch}/demultiplexed",
        otutab = lambda wildcards: get_profiles(wildcards, classifier = wildcards.classifier),
    output: "batches/{classifier}_{batch}.tsv"
    run:        
        import pandas as pd
        # merge tsv files
        otu_table = pd.DataFrame()
        for f in input.otutab:
            # get file name
            barcode = os.path.basename(f).split(".")[0]
            # read tsv file, first line is not header
            table = pd.read_csv(f, sep="\t", header = None)
            # rename columns
            table = table.rename(columns={0: "tax_id", 1: "taxonomy", 2: "estimated counts"})
            # if tax_id is 1, taxonomy is "unassigned"
            table.loc[table["tax_id"] == 1, "taxonomy"] = "unassigned"
            # use integer for "estimated counts"
            table["estimated counts"] = table["estimated counts"].astype(int)
            # rename "estimated counts" to sample name
            table = table.rename(columns={"estimated counts": barcode})
            # merge table into otu_table
            if otu_table.empty:
                otu_table = table
            else:
                otu_table = pd.merge(otu_table, table, on=["tax_id", "taxonomy"], how="outer")
        # fill NaN with 0
        otu_table = otu_table.fillna(0)
        # write otu_table to file; integer without decimal
        otu_table.to_csv(output[0], sep="\t", index=False, float_format="%.0f")
