localrules: guppy, demux_check, collect_fastq

rule guppy:
    # need to bind INPUT_DIR if not in workdir
    output:
        fastq_pass = temp(directory("{batch}/fastq_pass")),
        demux = temp(directory("{batch}/demux_guppy")),
    singularity: "docker://genomicpariscentre/guppy:3.3.3"
    log: "logs/demultiplex/guppy_{batch}.log"
    benchmark: "benchmarks/demultiplex/guppy_{batch}.txt"
    threads: config["threads"]["large"]
    params:
        fq = config["basecall_fq"],
        barcode_kits = config["guppy"]["barcode_kits"],
    resources:
        mem = config["mem"]["large"],
    shell: 
        """
        if [ {params.fq} == "None" ]; then
            echo "ERROR: 'basecall_fq' not found in config."
            exit 1
        fi
        mkdir -p {output.fastq_pass}
        cp {params.fq} {output.fastq_pass} -p
        guppy_barcoder -i {output.fastq_pass} -s {output.demux} -t {threads} --barcode_kits {params.barcode_kits} --trim_barcodes 2>{log}
        """
    
rule minibar:
    output:
        demux = temp(directory("{batch}/demux_minibar")),
    conda: "../envs/minibar.yaml"
    params:
        fq = config["basecall_fq"],
        args = config["minibar"]["args"],
    log: "logs/demultiplex/minibar_{batch}.log"
    benchmark: "benchmarks/demultiplex/minibar_{batch}.txt"
    resources:
        mem = config["mem"]["normal"],
        time = config["runtime"]["default"],
    shell: 
        """
        mkdir {output.demux} -p
        python {workflow.basedir}/scripts/minibar.py {workflow.basedir}/resources/data/index.txt {params.fq} -F -T -P {output.demux}/ {params.args} 2> {log}

        mkdir -p {output}/unclassified {output}/mult
        # create barcode list
        cut -f1 {workflow.basedir}/resources/data/index.txt | sed 1d > {output}/barcodes.txt 2>> {log}
        # sample in barcode list in a dir with batchid
        while read p; do
            if [ -f {output}/$p.fastq ]; then
                # if dir not exists, mkdir and move file
                if [ ! -d {output}/$p ]; then
                    mkdir {output}/$p 
                fi
                mv {output}/$p.fastq {output}/$p
            fi
        done < {output}/barcodes.txt
        # if file exists, mv to unclassified
        if [ -f {output}/unk.fastq ]; then
           mv {output}/unk.fastq {output}/unclassified
        fi
        # if .fastq file exists, mv to mult
        if [ -n "$(ls -A {output}/*.fastq 2>/dev/null)" ]; then
           mv {output}/*.fastq {output}/mult
        fi
        """

def check_demux_dir(dir_path = config["demultiplex_dir"]):
    if dir_path is not None:
        if not os.path.isdir(dir_path):
            raise ValueError("\n  'demultiplex_dir' ({}) in config not found.\n".format(dir_path))
        else:
            if not os.listdir(dir_path):
                raise ValueError("\n  'demultiplex_dir' ({}) in config is empty.\n".format(dir_path))
            else:
                # shall contain at least one barcode sub-directory, {barcode, [a-zA-Z]+[0-9]+}/{runid}.fastq
                glob_pattern = os.path.join(dir_path, "[a-zA-Z]*[0-9]*", "*.fastq")
                import glob
                if not glob.glob(glob_pattern):
                    raise ValueError("\n  'demultiplex_dir' ({}) in config does not contain barcode sub-directories.\n  A barcode folder contains unzipped fastq files, named as [a-zA-Z]+[0-9]+.\n".format(dir_path))

rule get_demux_external:
    output: temp(directory("{batch}/demux_external"))
    params:
        indir = config["demultiplex_dir"] 
    shell: "cp -r {params.indir} {output}"

def check_basecall_fq(fq = config["basecall_fq"]):
    if fq is not None:
        if not os.path.isfile(fq):
            raise ValueError("\n  'basecall_fq' ({}) in config not found.\n".format(fq))
        else:
            if os.path.getsize(fq) == 0:
                raise ValueError("\n  'basecall_fq' ({}) in config is empty.\n".format(fq))

# get demux input
def get_demux(demux=config["demuxer"], demux_external=config["demultiplex_dir"], batch_id=BATCH_ID):
    if demux_external is not None:
        check_demux_dir(demux_external)
        return batch_id + "/demux_external"
    else:
        if demux != "guppy" and demux != "minibar":
            raise ValueError("Demultiplexer not recognized. Choose guppy or minibar in config.")
        check_basecall_fq()
        return batch_id + "/demux_" + demux

checkpoint demux_check:
    input: get_demux()
    output: temp(directory("{batch}/demultiplexed"))
    log: "logs/demultiplex/check_{batch}.log"
    params:
        nreads_m=config["nreads_m"],
    shell: 
        """
        mv {input} {output} > {log} 2>&1
        # rm shallow sequencing to avoid bardcode bleeding
        mkdir {output}/suspected -p >> {log} 2>&1
        for i in {output}/*/
        do
            if [ "$i" == "{output}/suspected/" ] || [ "$i" == "{output}/unclassified/" ] || [ "$i" == "{output}/mult/" ]
            then
                continue
            fi
            nlines=$(cat $i/* | wc -l)
            nreads=$((nlines / 4))
            if [ $nreads -lt {params.nreads_m} ]
            then
                mv "$i" {output}/suspected/ >> {log} 2>&1
                echo "$i moved to {output}/suspected/ due to shallow sequencing < {params.nreads_m}" >> {log} 2>&1
            fi
        done
        """

rule collect_fastq:
    input:  "{batch}/demultiplexed/{barcode}"
    output: temp("{batch}/qc/{barcode}.fastq")
    log: "logs/demultiplex/collect_fastq/{barcode}_{batch}.log"
    benchmark: "benchmarks/demultiplex/collect_fastq/{barcode}_{batch}.txt"
    shell: "cat {input}/*.fastq | sed 's/^+.*/+/' | seqkit sana 2> {log} | seqkit rename -w0 -o {output} 2>> {log}"

def get_demux_barcodes(wildcards):
    barcodes = glob_wildcards(checkpoints.demux_check.get(**wildcards).output[0]
     + "/{barcode, [a-zA-Z]+[0-9]+}/{runid}.fastq").barcode
    return sorted(set(barcodes))
