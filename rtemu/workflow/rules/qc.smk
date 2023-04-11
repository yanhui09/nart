# primer info
fprimers = config["fprimer"]
rprimers = config["rprimer"]

def revcomp(seq):
    return seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))[::-1]
# nanopore possibly sequences either strand
def seqs_join(primer1, primer2):
    joined = '-g ' + primer1 + '...' + revcomp(primer2)
    return joined
def linked_pattern(primers1, primers2):
    primers1_values = list(primers1.values())
    primers2_values = list(primers2.values())
    linked = [seqs_join(primer1, primer2) for primer1 in primers1_values for primer2 in primers2_values]
    return ' '.join(linked)
# pattern
f5_pattern1 = linked_pattern(fprimers, rprimers)
f5_pattern2 = linked_pattern(rprimers, fprimers)

rule subsample:
    input: rules.collect_fastq.output
    output:
        p = temp("{batch}/qc/subsampled/{barcode}_p.fastq"),
        n = temp("{batch}/qc/subsampled/{barcode}.fastq"),
    conda: "../envs/seqkit.yaml"
    params:
        n = config["seqkit"]["n"],
    log: "logs/qc/subsample/{barcode}_{batch}.log"
    benchmark: "benchmarks/qc/subsample/{barcode}_{batch}.txt"
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["simple"],
    shell:
        """
        nlines=$(cat {input} | wc -l)
        nreads=$((nlines / 4))
        p=$(echo "scale=1; {params.n} / $nreads + 0.1" | bc)
        if (( $(echo "$p > 1" | bc -l) )); then
            p=1
        fi
        seqkit sample -p $p -j {threads} {input} -o {output.p} -w0 -s123 2> {log}
        seqkit head -n {params.n} -j {threads} {output.p} -o {output.n} -w0 2>> {log}
        """

def get_raw(subsample = config["subsample"], n = config["seqkit"]["n"]):
    check_val("subsample", subsample, bool)
    check_val("n[seqkit]", n, int)
    if subsample is True:
        return rules.subsample.output.n
    else:
        return rules.collect_fastq.output

# trim primers, process two strands differently
rule trim_primers:
    input: get_raw()
    output: 
        trimmed = temp("{batch}/qc/primers_trimmed/{barcode}F.fastq"),
        untrimmed = temp("{batch}/qc/primers_untrimmed/{barcode}F.fastq"),
    params:
        f = f5_pattern1,
        e = config["cutadapt"]["max_errors"],
        O = config["cutadapt"]["min_overlap"],
        m = 1,
    log: "logs/qc/trim_primersF/{barcode}_{batch}.log"
    benchmark: "benchmarks/qc/trim_primersF/{barcode}_{batch}.txt"
    threads: config["threads"]["large"]
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["simple"],
    shell:
        """
        cutadapt \
        -j {threads} \
        -e {params.e} -O {params.O} -m {params.m} \
        {params.f} \
        --untrimmed-output {output.untrimmed} \
        -o {output.trimmed} \
        {input} \
        > {log} 2>&1
        """

use rule trim_primers as trim_primersR with:
    input: 
        rules.trim_primers.output.untrimmed
    output:
        trimmed = temp("{batch}/qc/primers_trimmed/{barcode}R.fastq"),
        untrimmed = temp("{batch}/qc/primers_untrimmed/{barcode}.fastq"),
    params:
        f = f5_pattern2,
        e = config["cutadapt"]["max_errors"],
        O = config["cutadapt"]["min_overlap"],
        m = 1,
    log: 
        "logs/qc/trim_primersR/{barcode}_{batch}.log"
    benchmark: 
        "benchmarks/qc/trim_primersR/{barcode}_{batch}.txt"

# reverse complement for reverse strand
rule revcomp_fq:
    input: rules.trim_primersR.output.trimmed
    output: temp("{batch}/qc/primers_trimmed/{barcode}R_revcomp.fastq")
    log: "logs/qc/revcomp_fq/{barcode}_{batch}.log"
    benchmark: "benchmarks/qc/revcomp_fq/{barcode}_{batch}.txt"
    threads: config["threads"]["normal"]
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["simple"],
    shell: "seqkit seq -j {threads} -r -p -t dna {input} > {output} 2> {log}"

# option to trim or not
def trim_check(trim = config["trim"], subsample = config["subsample"], n = config["seqkit"]["n"]):
    check_val("trim", trim, bool)
    out = [rules.trim_primers.output.trimmed, rules.revcomp_fq.output]
    if trim is False:
        out = get_raw(subsample, n)
    return out

rule q_filter:
    input: trim_check()
    output: temp("{batch}/qc/qfilt/{barcode}.fastq")
    params:
        Q = config["seqkit"]["min_qual"],
        m = config["seqkit"]["min_len"],
        M = config["seqkit"]["max_len"],
    log: "logs/qc/q_filter/{barcode}_{batch}.log"
    benchmark: "benchmarks/qc/q_filter/{barcode}_{batch}.txt"
    threads: config["threads"]["normal"]
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["simple"],
    shell: "cat {input} | seqkit seq -j {threads} -Q {params.Q} -m {params.m} -M {params.M} -i > {output} 2> {log}"

checkpoint exclude_empty_fqs:
    input: lambda wc: expand("{batch}/qc/qfilt/{barcode}.fastq", barcode=get_demux_barcodes(wc))
    output: temp(directory("{batch}/qc/qfilt/empty"))
    run:
        import shutil
        if not os.path.exists(output[0]):
            os.makedirs(output[0])
        for i in list(input):
            if os.stat(i).st_size == 0:
                shutil.move(i, output[0])

def get_qced_barcodes(wildcards, batch_id = BATCH_ID):
    barcodes = get_demux_barcodes(wildcards)
    barcodes_empty = glob_wildcards(checkpoints.exclude_empty_fqs.get(**wildcards).output[0] + "/{barcode}.fastq").barcode
    barcodes_empty = sorted(set(barcodes_empty))
    barcodes = [b for b in barcodes if b not in barcodes_empty]
    return barcodes

def get_filt(wildcards):
    barcodes = get_qced_barcodes(wildcards) 
    return expand("{{batch}}/qc/qfilt/{barcode}.fastq", barcode=barcodes)