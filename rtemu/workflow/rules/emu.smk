rule emu:
    input: 
        rules.emu_prebuilt.output,
        fq = "{batch}/qc/qfilt/{barcode}.fastq",
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

def get_emu(wildcards):
    barcodes = get_qced_barcodes(wildcards) 
    return expand("{{batch}}/emu/{barcode}_rel-abundance.tsv", barcode=barcodes)

localrules: emu_merge
rule emu_merge:
    input: 
        demux_dir = "{batch}/demultiplexed",
        otutab = get_emu,
    output: "otu_table_{batch}.tsv"
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
            # combine "phylum", "class", "order", "family", "genus", "species" into "taxonomy"
            table["taxonomy"] = table["phylum"] + ";" + table["class"] + ";" + table["order"] + ";" + table["family"] + ";" + table["genus"] + ";" + table["species"]
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







