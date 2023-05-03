from email.policy import default
import click
import os
import subprocess
from .log import logger
from .config import init_conf

from snakemake import load_configfile
from . import __version__

def handle_max_mem(max_mem, profile):
    "Specify maximum memory (GB) to use. Memory is controlled by profile in cluster execution."
    "For numbers >1 its the memory in GB. "
    "For numbers <1 it's the fraction of available memory."

    if profile is not None:

        if max_mem is not None:
            logger.info(
                "Memory requirements are handled by the profile, I ignore max-mem argument."
            )
        # memory is handled via the profile, user should know what he is doing
        return ""
    else:
        import psutil
        from math import floor

        # calulate max  system meory in GB (float!)
        max_system_memory = psutil.virtual_memory().total / (1024**3)

        if max_mem is None:
            max_mem = 0.95
        if max_mem > 1:

            if max_mem > max_system_memory:
                logger.critical(
                    f"You specified {max_mem} GB as maximum memory, but your system only has {floor(max_system_memory)} GB"
                )
                sys.exit(1)

        else:

            max_mem = max_mem * max_system_memory

        # specify max_mem_string including java mem and max mem

        return f" --resources mem={floor(max_mem)} mem_mb={floor(max_mem*1024)} java_mem={floor(0.85* max_mem)} "

def get_snakefile(file="workflow/Snakefile"):
    sf = os.path.join(os.path.dirname(os.path.abspath(__file__)), file)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file; tried %s" % sf)
    return sf

def run_smk(workflow, workdir, configfile, jobs, maxmem, profile, dryrun, snake_args, snakefile, exit_on_error, suppress):
    """
    Start NAWF in a single batch.
    Most snakemake arguments can be appended, for more info see 'snakemake --help'
    """
    if not suppress:
        logger.info(f"NAWF version: {__version__}")
    if not os.path.exists(configfile):
        logger.critical(f"Workflow config file not found: {configfile}\nGenerate a config file using 'nawf config'")
        exit(1)
    
    conf = load_configfile(configfile)
    # if snake_args (tuple, whitespace removed) contain "basecalled_dir=", take everything after as basecalled_dir
    for arg in snake_args:
        if "basecall_fq=" in arg:
            conf["basecall_fq"] = arg.split("=")[1]
    basecall_fq = conf["basecall_fq"]
    db_dir = conf["database_dir"]
    cmd = (
        "snakemake "
        "{wf} "
        "--directory '{workdir}' "
        "--snakefile '{snakefile}' "
        "--configfile '{configfile}' "
        "--use-conda {conda_prefix} "
        "{singularity_prefix} "
        "{singularity_args} "
        "{dryrun} "
        "--rerun-triggers mtime --rerun-incomplete --scheduler greedy "
        "--jobs {jobs} --nolock "
        " {max_mem} "
        " {profile} "
        " {args} "
    ).format(
        wf=workflow if workflow is not None else "",
        workdir=workdir,
        snakefile=snakefile,
        configfile=configfile,
        conda_prefix="--conda-prefix '" + os.path.join(db_dir, "conda_envs") + "'",
        singularity_prefix="--use-singularity --singularity-prefix '" + os.path.join(db_dir, "singularity_envs") + "'"
        if basecall_fq is not None else "",
        singularity_args="--singularity-args '--bind " +
        os.path.dirname(snakefile) + "/resources/guppy_barcoding/:/opt/ont/guppy/data/barcoding/," + basecall_fq + "'"
        if basecall_fq is not None else "",
        dryrun="--dryrun" if dryrun else "",
        jobs=int(jobs) if jobs is not None else 1,
        max_mem=handle_max_mem(maxmem, profile),
        profile="" if (profile is None) else "--profile {}".format(profile),
        args=" ".join(snake_args),
    )
    if not suppress:
        logger.debug("Executing: %s" % cmd)
    
    try:
        if suppress:
            subprocess.check_call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        else:
            subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        logger.critical(e)
        if exit_on_error is True:
            exit(1)

# custom alo(at least one) class with mutex (mutually exclusive) classs
# https://stackoverflow.com/questions/44247099/click-command-line-interfaces-make-options-required-if-other-optional-option-is/
class AloMutex(click.Option):
    def __init__(self, *args, **kwargs):
        self.required_if_not:list = kwargs.pop("required_if_not")
        self.not_required_if:list = kwargs.pop("not_required_if")

        if self.required_if_not:
            kwargs["help"] = (
                kwargs.get("help", "") + " Option is required if '" + "', '".join(self.required_if_not) + "' not provided."
                ).strip()
        if self.not_required_if:
            kwargs["help"] = (
                kwargs.get("help", "") + " Option is mutually exclusive with '" + "', '".join(self.not_required_if) + "'."
                ).strip()
        super(AloMutex, self).__init__(*args, **kwargs)

    def handle_parse_result(self, ctx, opts, args):
        current_opt:bool = self.name in opts
        if self.required_if_not and all(alo_opt not in opts for alo_opt in self.required_if_not) and not current_opt:
            raise click.UsageError(
                "at least one of '" + "', '".join(self.required_if_not) + 
                "' and '" + str(self.name) +  "' options is provided."
                )
        if self.not_required_if:
            for mutex_opt in self.not_required_if:
                if mutex_opt in opts:
                    if current_opt:
                        raise click.UsageError(
                            "'" + str(self.name) + "' is mutually exclusive with '" + str(mutex_opt) + "'."
                            )
                    else:
                        self.prompt = None
        return super(AloMutex, self).handle_parse_result(ctx, opts, args)

@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(
    __version__,
    "-v",
    "--version",
    )
@click.pass_context
def cli(self):
    """
    NAWF: A sub-tool to run Nanopore Amplicon WorkFlow.
    The workflow command initiates the NAWF in a single batch,
    using either a fastq file from one ONT run or a fastq file generated during sequencing.
    To follow updates and report issues, see: https://github.com/yanhui09/nart.
    """
    pass

@cli.command(
    'run',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Start workflow in a single batch.'
)
@click.option(
    "-w",
    "--workdir",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="Workflow working directory.",
    show_default=True,
    default=".",
)
@click.option(
    "-c",
    "--configfile",
    type=click.Path(exists=True, dir_okay=False, resolve_path=True),
    help="Workflow config file. Use config.yaml in working directory if not specified.",
    default=None,
)
@click.option(
    "-j",
    "--jobs",
    type=int,
    default=6,
    show_default=True,
    help="Maximum jobs to run in parallel.",
)
@click.option(
    "-m",
    "--maxmem",
    type=float,
    default=None,
    show_default=True,
    help=handle_max_mem.__doc__,
)
@click.option(
    "--profile",
    default=None,
    help="Snakemake profile for cluster execution.",
)
@click.option(
    "-n",
    "--dryrun",
    is_flag=True,
    default=False,
    show_default=True,
    help="Dry run.",
)

@click.argument(
    "workflow",
    type=click.Choice(
        ["init","demux", "qc","all"]
    ),
)
@click.argument("snake_args", nargs=-1, type=click.UNPROCESSED)
def run_workflow(workflow, workdir, configfile, jobs, maxmem, profile, dryrun, snake_args):
    """
    Run NAWF in a single batch.
    """
    sf = "workflow/Snakefile" 
    snakefile = get_snakefile(sf)
    configfile_run = os.path.join(workdir, "config.yaml") if configfile is None else configfile
    run_smk(workflow, workdir, configfile_run, jobs, maxmem, profile, dryrun, snake_args, snakefile, exit_on_error=True, suppress=False)

# workflow config
# initialize config file
@cli.command(
    'config',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Generate the workflow config file.',
)
@click.option(
    "-b",
    "--bascfq",
    help="Path to a basecalled fastq file.",
    type=click.Path(dir_okay=True,writable=True,resolve_path=True),
    cls = AloMutex,
    required_if_not = [],
    not_required_if = ["demuxdir"],
)
@click.option(
    "-x",
    "--demuxdir",
    help="Path to a directory of demultiplexed fastq files.",
    type=click.Path(dir_okay=True,writable=True,resolve_path=True),
    cls = AloMutex,
    required_if_not = [],
    not_required_if = ["bascfq"],
)
@click.option(
    "-d",
    "--dbdir",
    help="Path to the taxonomy databases.",
    type=click.Path(dir_okay=True,writable=True,resolve_path=True),
    required=True
)
@click.option(
    "-w",
    "--workdir",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="Output directory for NAWF.",
    show_default=True,
    default=".",
)
@click.option(
    "--demuxer",
    type=click.Choice(["guppy", "minibar"]),
    default="guppy",
    show_default=True,
    help="Demultiplexer.",
)
@click.option(
    "--fqs-min",
    type=int,
    default=50,
    show_default=True,
    help="Minimum number of reads for the demultiplexed fastqs.",
)
@click.option(
    "--subsample",
    is_flag=True,
    default=False,
    show_default=True,
    help="Subsample the reads.",
)
@click.option(
    "--trim",
    is_flag=True,
    default=False,
    show_default=True,
    help="Trim primers.",
)
@click.option(
    "--classifier",
    type=click.Choice(["emu", "minimap2lca", "blast2lca"]),
    default="emu",
    show_default=True,
    help="Classifier.",
)
@click.option(
    "--jobs-min",
    type=int,
    default=2,
    show_default=True,
    help="Number of jobs for common tasks.",
)
@click.option(
    "--jobs-max",
    type=int,
    default=6,
    show_default=True,
    help="Number of jobs for threads-dependent tasks.",
)
def config_workflow(
    bascfq, demuxdir, dbdir, workdir, demuxer, fqs_min, subsample, trim, 
    classifier, jobs_min, jobs_max):
    """
    Config NAWF.
    """ 
    logger.info(f"NAWF version: {__version__}")
    init_conf(
        bascfq, demuxdir, dbdir, workdir, "config.yaml", demuxer, fqs_min, subsample,
        trim, classifier, jobs_min, jobs_max)
   
if __name__ == "__main__":
    cli()
