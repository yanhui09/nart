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

# filter chimeras with yacrd
rule minimap2ava:
    input: get_raw()
    output: temp("{batch}/qc/yacrd/{barcode}.paf")
    conda: "../envs/yacrd.yaml"
    params:
        x = "ava-ont",
        g = 500,
        f = 1000,
    log: "logs/qc/yacrd/{barcode}_{batch}_ava.log"
    benchmark: "benchmarks/qc/yacrd/{barcode}_{batch}_ava.txt"
    threads: config["threads"]["large"]
    resources:
        mem = config["mem"]["large"],
        time = config["runtime"]["simple"],
    shell: "minimap2 -x {params.x} -g {params.g} -f {params.f} -t {threads} {input} {input} > {output} 2> {log}"

rule yacrd:
    input: 
        fq = get_raw(),
        ava = rules.minimap2ava.output
    output: temp("{batch}/qc/yacrd/{barcode}.fastq")
    conda: "../envs/yacrd.yaml"
    params:
        c = config["yacrd"]["c"],
        n = config["yacrd"]["n"],
    log: "logs/qc/yacrd/{barcode}_{batch}_scrubb.log"
    benchmark: "benchmarks/qc/yacrd/{barcode}_{batch}_scrubb.txt"
    threads: config["threads"]["large"]
    resources:
        mem = config["mem"]["large"],
        time = config["runtime"]["simple"],
    shell: "yacrd -i {input.ava} -o {log} -c {params.c} -n {params.n} -t {threads} filter -i {input.fq} -o {output} 2>> {log}"

def get_chimera_free(chimera_filt= config["chimera_filt"]):
    check_val("chimera_filt", chimera_filt, bool)
    if chimera_filt is True:
        return rules.yacrd.output
    else:
        return get_raw()

# check primer-pattern, process two strands independently
rule check_primers:
    input: get_chimera_free()
    output: 
        passed = temp("{batch}/qc/primers_passed/{barcode}F.fastq"),
        unpassed = temp("{batch}/qc/primers_unpassed/{barcode}F.fastq"),
    params:
        f = f5_pattern1,
        e = config["cutadapt"]["max_errors"],
        O = config["cutadapt"]["min_overlap"],
        m = 1,
        action = config["cutadapt"]["action"],
    log: "logs/qc/check_primersF/{barcode}_{batch}.log"
    benchmark: "benchmarks/qc/check_primersF/{barcode}_{batch}.txt"
    threads: config["threads"]["large"]
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["simple"],
    shell:
        """
        cutadapt \
        --action={params.action} \
        -j {threads} \
        -e {params.e} -O {params.O} -m {params.m} \
        {params.f} \
        --untrimmed-output {output.unpassed} \
        -o {output.passed} \
        {input} \
        > {log} 2>&1
        """

use rule check_primers as check_primersR with:
    input: 
        rules.check_primers.output.unpassed
    output:
        passed = temp("{batch}/qc/primers_passed/{barcode}R.fastq"),
        unpassed = temp("{batch}/qc/primers_unpassed/{barcode}.fastq"),
    params:
        f = f5_pattern2,
        e = config["cutadapt"]["max_errors"],
        O = config["cutadapt"]["min_overlap"],
        m = 1,
        action = config["cutadapt"]["action"],
    log: 
        "logs/qc/check_primersR/{barcode}_{batch}.log"
    benchmark: 
        "benchmarks/qc/check_primersR/{barcode}_{batch}.txt"

# reverse complement for reverse strand
rule revcomp_fq:
    input: rules.check_primersR.output.passed
    output: temp("{batch}/qc/primers_passed/{barcode}R_revcomp.fastq")
    log: "logs/qc/revcomp_fq/{barcode}_{batch}.log"
    benchmark: "benchmarks/qc/revcomp_fq/{barcode}_{batch}.txt"
    threads: config["threads"]["normal"]
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["simple"],
    shell: "seqkit seq -j {threads} -r -p -t dna {input} > {output} 2> {log}"

def primer_check(primer_check = config["primer_check"], subsample = config["subsample"], n = config["seqkit"]["n"]):
    check_val("primer_check", primer_check, bool)
    out = [rules.check_primers.output.passed, rules.revcomp_fq.output]
    if primer_check is False:
        out = get_raw(subsample, n)
    return out

rule q_filter:
    input: primer_check()
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