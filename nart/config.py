from .log import logger

import os
from ruamel.yaml import YAML

def init_conf(
    bascfq,
    demuxdir,
    dbdir,
    workdir,
    config="config.yaml",
    demuxer = "guppy",
    nreads_m = 1000,
    subsample=False,
    chimera_filt=False,
    primer_check=False,
    classifier="minimap2lca",
    jobs_m=2,
    jobs_M=6,
):
    """
    Reads template config file with comments from ./template_config.yaml
    updates it by the parameters provided.
    Args:
        bascfq (str): path to a directory of basecalled fastq files
        demuxdir (str): path to a directory of demultiplexed fastq files
        dbdir (str): path to the taxonomy database
        workdir (str): path to the working directory
        config (str): the config filename
        demuxer (str): the demultiplexer [default: "guppy"]
        nreads_m (int): minimum number of reads for the demultiplexed fastqs
        subsample (bool): if True, subsample the reads [default: False]
        chimera_filt (bool): if False, do not filter possible chimeric reads [default: False]
        primer_check (bool): if False, do not check primer pattern [default: False]
        minimap2lca (str): the classifier [default: "minimap2lca"]
        jobs_m (int): number of jobs for common tasks [default: 2]
        jobs_M (int): number of jobs for threads-dependent tasks [default: 6]
   """
    os.makedirs(dbdir, exist_ok=True)
    os.makedirs(workdir, exist_ok=True)
    
    yaml = YAML()
    template_conf_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "workflow/config_template.yaml")
    
    with open(template_conf_file) as template_conf:
        conf = yaml.load(template_conf)
    
    conf["basecall_fq"] = bascfq
    conf["demultiplex_dir"] = demuxdir    
    conf["database_dir"] = dbdir
    conf["demuxer"] = demuxer
    conf["nreads_m"] = nreads_m
    conf["subsample"] = subsample
    conf["chimera_filt"] = chimera_filt
    conf["primer_check"] = primer_check
    conf["classifier"] = classifier
    conf["threads"]["normal"] = jobs_m
    conf["threads"]["large"] = jobs_M
    
    if os.path.exists(os.path.join(workdir, config)):
        logger.warning(f"Config file [{config}] already exists in {workdir}.")
    else:
        with open(os.path.join(workdir, config), "w") as conf_file:
            yaml.dump(conf, conf_file)
        logger.info(f"Config file [{config}] created in {workdir}.")
    