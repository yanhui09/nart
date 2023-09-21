def get_emu_database(spikein=config["spikein_fasta"]):
    if spikein == "none":
        return rules.emu_prebuilt.output
    else:
        return rules.emu_spikein.output

rule emu:
    input: 
        get_emu_database(),
        fq = rules.q_filter.output,
    output: temp("{batch}/emu/{barcode}_rel-abundance.tsv")
    conda: "../envs/emu.yaml"
    params: 
        db = DATABASE_DIR + "/emu/" + DATABASE_PREBUILT + "_prebuilt" if config["spikein_fasta"] == "none" else DATABASE_DIR + "/emu/" + DATABASE_PREBUILT + "_prebuilt_spikein",
        outdir = os.getcwd() + "/{batch}/emu",
        min_abundance = config["emu"]["min_rel_abundance"],
    log: "logs/emu/{barcode}_{batch}.log"
    benchmark: "benchmarks/emu/{barcode}_{batch}.txt"
    threads: config["threads"]["large"]
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["default"],
    shell:
        """
        emu abundance --db {params.db} {input.fq} --min-abundance {params.min_abundance} --keep-counts --threads {threads} --output-dir {params.outdir} > {log} 2>&1 || touch {output}
        # if {output} empty, rm *_emu_alignments.sam and send warning message to log
        if [ ! -s {output} ]
        then
            rm {params.outdir}/{wildcards.barcode}_emu_alignments.sam -f
            echo "EMU failed for {wildcards.barcode}, pseduo abundacne file was created." >> {log}
        fi
        """

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
        threshold = config["emu"]["min_rel_abundance"],
        export_rel_abundance = config["emu"]["export_rel_abundance"],
    resources:
        mem = config["mem"]["normal"],
    run:
        import pandas as pd
        def merge_table(list_otutab, col_abundance, file_out):
            # merge tsv files
            otu_table = pd.DataFrame()
            for f in list_otutab:
                if os.stat(f).st_size == 0:
                    continue
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
                # if taxonomy is ";;;;;;spikein", taxonomy is "spikein"
                table.loc[table["taxonomy"] == ";;;;;;spikein", "taxonomy"] = "spikein"
                # only keep "tax_id","taxonomy", type_abundance = "estimated counts" or "abundance"
                table = table[["tax_id", "taxonomy", col_abundance]]
                if col_abundance == "estimated counts":
                    # use integer for "estimated counts"
                    table[[col_abundance]] = table[[col_abundance]].astype(int)
                else:
                    # use float for "abundance"
                    table[[col_abundance]] = table[[col_abundance]].astype(float)
                # rename col_abundance to sample name
                table = table.rename(columns={col_abundance: barcode})
                # merge table into otu_table
                if otu_table.empty:
                    otu_table = table
                else:
                    otu_table = pd.merge(otu_table, table, on=["tax_id", "taxonomy"], how="outer")
            # fill NaN with 0
            otu_table = otu_table.fillna(0)
            # if dataframe is empty, use open() to create empty file to avoid an empty line in output
            if otu_table.empty:
                open(file_out, 'w').close()
            else:
                # write otu_table to file; integer without decimal
                if col_abundance == "estimated counts":
                    otu_table.to_csv(file_out, sep="\t", index=False, float_format="%.0f")
                else:
                    otu_table.to_csv(file_out, sep="\t", index=False)
            
        # use estimated counts as default
        merge_table(input.otutab, "estimated counts", output[0])
        if params.export_rel_abundance is True:
            emu_extra = os.getcwd() + "/emu_extra/"
            os.makedirs(emu_extra, exist_ok=True)
            merge_table(input.otutab, "abundance", emu_extra + wildcards.batch + "_rel.tsv")
        # possibly excluded
        otutab_unpassed = [os.path.dirname(i) + "/" + os.path.basename(i).replace("-abundance", "-abundance-threshold-" + str(params.threshold)) for i in input.otutab]
        # if params.export_rel_abundance is True or any(otutab_unpassed) exists, create emu_extra folder
        if any(os.path.exists(i) for i in otutab_unpassed):
            emu_extra = os.getcwd() + "/emu_extra/"
            os.makedirs(emu_extra, exist_ok=True)
            merge_table(otutab_unpassed, "estimated counts", emu_extra + wildcards.batch + "_excluded.tsv")
            if params.export_rel_abundance is True:
                merge_table(otutab_unpassed, "abundance", emu_extra + wildcards.batch + "_excluded_rel.tsv")

def get_silva_database(mode="minimap2", spikein=config["spikein_fasta"], taxmap=False):
    if spikein == "none":
        silva = "silva"
    else:
        silva = "silva_spikein"
    if taxmap is True:
        return DATABASE_DIR + "/{silva}/silva_to_ncbi.map".format(silva=silva)
    else:
        if mode == "minimap2":
            return DATABASE_DIR + "/{silva}/SILVA_ssu_nr99_filt.mmi".format(silva=silva)
        elif mode == "blastn":
            return  multiext(DATABASE_DIR + "/{silva}/SILVA_ssu_nr99_filt.fasta".format(silva=silva),
                ".ndb",
                ".nhr",
                ".nin",
                ".not",
                ".nsq",
                ".ntf",
                ".nto"
                )
        else:
            raise ValueError("mode must be minimap2 or blastn")

# minimap2lca
rule minimap2:
    input:
        mmi = get_silva_database(mode="minimap2"),
        fq = rules.q_filter.output,
    output: 
        sam4 = temp("{batch}/minimap2lca/{barcode}.sam4"),
        sam = temp("{batch}/minimap2lca/{barcode}.sam")
    conda: "../envs/lca.yaml"
    params:
        K = '500M',
        x = 'map-ont',
        f = 10000,
        # https://lh3.github.io/minimap2/minimap2.html#10
        max_de = 0.1,
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["default"],
    log: "logs/minimap2lca/{barcode}_{batch}_minimap2.log"
    benchmark: "benchmarks/minimap2lca/{barcode}_{batch}_minimap2.txt"
    threads: config["threads"]["large"]
    shell:
        """
        if [ ! -s {input.fq} ]; then
            touch {output.sam4} {output.sam}
        else
            minimap2 -t {threads} -K {params.K} -ax {params.x} -f {params.f} --secondary=no {input.mmi} {input.fq} > {output.sam4} 2> {log}
            grep "^@" {output.sam4} > {output.sam}
            grep -v "^@" {output.sam4} | awk -F '\t|de:f:' '$(NF-1) < {params.max_de}' >> {output.sam} 2>> {log}
        fi
        """

rule minimap2rma:
    input:
        fq = rules.q_filter.output,
        sam = rules.minimap2.output.sam,
        silva2ncbi_map = get_silva_database(taxmap=True),
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
        if [ ! -s {input.fq} ] && [ ! -s {input.sam} ]; then
            touch {output}
        else
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
        fi
        """

rule blastn:
    input:
        get_silva_database(mode="blastn"), 
        db = get_silva_database(mode="blastn")[0].split(".ndb")[0],
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
        if [ ! -s {input.fq} ]; then
            touch {output.fa} {output.tab}
        else
            seqkit fq2fa {input.fq} -o {output.fa} > {log} 2>&1
            blastn -query {output.fa} -db {input.db} -out {output.tab} -outfmt 6 -num_threads {threads} -max_target_seqs {params.max_target_seqs} >> {log} 2>&1
        fi
        """

rule blast2rma:
    input:
        fa = rules.blastn.output.fa,
        tab = rules.blastn.output.tab,
        silva2ncbi_map = get_silva_database(taxmap=True),
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
        if [ ! -s {input.fa} ] && [ ! -s {input.tab} ]; then
            touch {output}
        else
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
        fi
        """

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
        if [ ! -s {input} ]; then
            touch {output.ncbi} {output.path} {output.ri}
        else
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
        fi
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
            if os.path.getsize(f) == 0:
                continue
            barcode = os.path.basename(f).split(".")[0]
            # read tsv file, first line is not header
            table = pd.read_csv(f, sep="\t", header = None)
            # rename columns
            table = table.rename(columns={0: "tax_id", 1: "taxonomy", 2: "estimated counts"})
            # if tax_id is 1, taxonomy is "unassigned"
            table.loc[table["tax_id"] == 1, "taxonomy"] = "unassigned"
            # if tax_id == 10710, taxonomy is "spikein"
            table.loc[table["tax_id"] == 10710, "taxonomy"] = "spikein"
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
        # if dataframe is empty, use open() to create empty file to avoid an empty line in output
        if otu_table.empty:
            open(output[0], 'w').close()
        else:
            # write otu_table to file; integer without decimal
            otu_table.to_csv(output[0], sep="\t", index=False, float_format="%.0f")
