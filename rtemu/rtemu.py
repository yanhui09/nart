import click
from .workflow import handle_max_mem, get_snakefile, run_smk
from . import __version__
import os
import time
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
import logging
logging.getLogger("watchdog.observers.inotify_buffer").setLevel(logging.WARNING)
import pandas as pd
from .app import run_server

def monitor_directory_for_new_files(directory_path, file_extension, timeout_seconds, output_file):
    """
    Monitors a directory for new files with a specific extension and writes the absolute path of each new file
    to a text file.

    :param directory_path: The path of the directory to monitor.
    :param file_extension: The file extension to monitor for (e.g. ".fastq").
    :param timeout_seconds: The number of seconds to wait before exiting if no new files are created.
    :param output_file: The path of the file to write the absolute path of each new file to.
    """
    class MyHandler(FileSystemEventHandler):
        def __init__(self, output_file):
            self.last_file_time = time.time()
            self.output_file = output_file

        def on_created(self, event):
            if event.is_directory:
                return None
            elif event.src_path.endswith(file_extension):
                print(f"New file created: {event.src_path}")
                with open(self.output_file, "a") as f:
                    f.write(os.path.abspath(event.src_path) + "\n")
                self.last_file_time = time.time()

        def on_modified(self, event):
            if event.is_directory:
                return None
            elif event.src_path.endswith(file_extension):
                self.last_file_time = time.time()

    event_handler = MyHandler(output_file)
    observer = Observer()
    observer.schedule(event_handler, path=directory_path, recursive=False)
    observer.start()

    try:
        while True:
            if time.time() - event_handler.last_file_time > timeout_seconds:
                print(f"No new files with extension {file_extension} created in the last {timeout_seconds} seconds. Exiting.")
                break
            time.sleep(1)
    except KeyboardInterrupt:
        observer.stop()
    observer.join()

# rtemu monitor
# minoitor the new fqs in a folder, update the fqs.txt and run the emuwf
# to do --visual render a html report
def run_monitor(query, extension, workdir, timeout_seconds, output_file="fqs.txt"):
    # check if the query is a directory
    if not os.path.isdir(query):
        raise ValueError("The query must be a directory.")
    # check workdir if exists
    if not os.path.exists(workdir):
        os.makedirs(workdir)
 
    output_path = workdir + "/" + output_file
    # if file exists, remove it
    if os.path.exists(output_path):
        os.remove(output_path)
    # create a new file 
    # and list the file with required extensions in the query folder if possible
    with open(output_path, "w") as f:
        for file in os.listdir(query):
            if file.endswith(extension):
                f.write(os.path.abspath(query + "/" + file) + "\n")                
    monitor_directory_for_new_files(query, extension, timeout_seconds, output_path)
 
@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(
    __version__,
    "-v",
    "--version",
    )
@click.pass_context
def cli(self):
    """
    RT-Emu: Real-time ONT amplicon workflow with Emu.
    To follow updates and report issues, see: https://github.com/yanhui09/rtemu.
    """
    pass

# one thread to monitor the new fqs
@cli.command(
    "monitor",
    context_settings=dict(ignore_unknown_options=True),
    short_help='Start RT-Emu to monitor a directory.'
)
@click.option(
    "-q",
    "--query",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="A query directory to monitor the new fastq files.",
    default=None,
)
@click.option(
    "-e",
    "--extension",
    type=str,
    default=".fastq",
    show_default=True,
    help="The file extension to monitor for (e.g. '.fastq').",
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
    "-t",
    "--timeout",
    type=int,
    default=30,
    show_default=True,
    help="Stop query if no new files were generated within the give minutes.",
)
def monitor(query, extension, workdir, timeout):
    """
    Start RT-Emu monitor.
    """
    run_monitor(query, extension, workdir, int(timeout)*60)
    exit(0)

# run
def merge_table(out_table, otu_table):
    # combine tables by "tax_id" and "taxonomy" if they are the same
    out_table = pd.merge(out_table, otu_table, on=["tax_id", "taxonomy"], how="outer", suffixes=("", "_x")).fillna(0)
    # add the values of "_x" to the original columns, and drop the "_x" columns
    col_x = [col for col in out_table.columns if col.endswith("_x")]
    for col in col_x:
        # without "_x"
        col_nox = col.replace("_x", "")
        out_table[col_nox] = out_table[col_nox] + out_table[col]
        out_table = out_table.drop(col, axis=1)
    return out_table    

def update_table(table_dir, out_table_path):
    ''''
    merge all otu tables in the table_dir as one;
    update table when out_table is newer than the otu batch table
    '''
    # list all .tsv files in the table_dir
    otu_tables = os.listdir(table_dir)
    otu_tables = [table_dir + "/" + table for table in otu_tables if table.endswith(".tsv")]
    # if the out_table exists, load it or create a new one
    if os.path.exists(out_table_path):
        out_table = pd.read_csv(out_table_path, sep="\t")
        for otu_table in otu_tables:
            # if the otu_table is newer than the out_table, update the out_table
            if os.path.getmtime(otu_table) > os.path.getmtime(out_table_path):
                otu_table = pd.read_csv(otu_table, sep="\t")
                out_table = merge_table(out_table, otu_table)
    else:
        # merge all otu tables
        out_table = pd.read_csv(otu_tables[0], sep="\t")
        for otu_table in otu_tables[1:]:
            otu_table = pd.read_csv(otu_table, sep="\t")
            out_table = merge_table(out_table, otu_table)
    # fill the nan with 0
    out_table = out_table.fillna(0)
    # write otu_table to file; integer without decimal
    out_table.to_csv(out_table_path, sep="\t", index=False, float_format="%.0f")
    
def run_rtemu(file_list, wait_minutes, workflow, workdir, configfile, jobs, maxmem, profile, dryrun, snake_args):
    # check if the workdir
    if not os.path.exists(workdir):
        raise ValueError("The workdir does not exist.")
    # if no file_list 
    if not os.path.exists(file_list):
        raise ValueError("The {} does not exist.".format(file_list))
    # stop when "fqs.txt" is not updated
    while True:
        # record file revised time
        file_time = os.path.getmtime(file_list)    
        # load the fqs.txt
        with open(file_list, "r") as f:
            file_paths = f.readlines()
            for file_path in file_paths:
                # remove the "\n" in the end of each line
                file_path = file_path.strip()
                # if the otu batch file exists, skip it
                batch_id = os.path.basename(file_path).split(".")[0]
                if os.path.exists(workdir + "/batches/" + batch_id + ".tsv"):
                    update_table(workdir + "/batches", workdir + "/otu_table.tsv")
                    continue
                snake_args_batch = ["--config", "basecall_fq=" + str(file_path)]
                snake_args = list(snake_args) + snake_args_batch
                snake_args = tuple(snake_args)
                sf = "workflow/Snakefile" 
                snakefile = get_snakefile(sf)
                configfile_run = os.path.join(workdir, "config.yaml") if configfile is None else configfile
                run_smk(workflow, workdir, configfile_run, jobs, maxmem, profile, dryrun, snake_args, snakefile, exit_on_error=True, suppress=False)
                update_table(workdir + "/batches", workdir + "/otu_table.tsv")
                if dryrun:
                    exit(0)
                
        # check if the fqs.txt is updated
        if file_time == os.path.getmtime(file_list):
            # if not updated, wait 30s and check again
            time.sleep(30)
            if file_time == os.path.getmtime(file_list):
                time.sleep(60*int(wait_minutes))
                if file_time == os.path.getmtime(file_list):
                    print("No new files added in the last {} minutes. Exiting.".format(wait_minutes))
                    break      
  
@cli.command(
    "run",
    context_settings=dict(ignore_unknown_options=True),
    short_help='Start RT-Emu workflow.'
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
    "-t",
    "--timeout",
    type=int,
    default=10,
    show_default=True,
    help="Stop run if no new files were updated in list within the given minutes.",
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
@click.argument("snake_args", nargs=-1, type=click.UNPROCESSED)
def run(workdir, timeout, configfile, jobs, maxmem, profile, dryrun, snake_args):
    """
    Start RT-Emu.
    """
    run_rtemu(workdir + "/fqs.txt", timeout, "all", workdir, configfile, jobs, maxmem, profile, dryrun, snake_args)
    exit(0)

# visual app
@cli.command(
    "visual",
    context_settings=dict(ignore_unknown_options=True),
    short_help='Start RT-Emu app to interactively visualize the results.'
)
@click.option(
    '-p', 
    '--port', 
    type=int,
    default=5000, 
    show_default=True,
    help='Port to run the app on.'
    )
@click.option(
    '-i', 
    '--input', 
    type=click.Path(dir_okay=True,writable=True,resolve_path=True),
    default='.', 
    show_default=True,
    help='Path to the working directory.'
    )
@click.option(
    '-w', 
    '--wait-time', 
    type=int,
    default=5, 
    show_default=True,
    help='Time to wait (in minutes) if input file is missing.'
    )
@click.option(
    "--relative",
    is_flag=True,
    default=False,
    show_default=True,
    help="Use relative abundance instead of absolute abundance.",
)
@click.option(
    "--rm-unmapped",
    is_flag=True,
    default=False,
    show_default=True,
    help="Remove unmapped reads from the table.",
)
@click.option(
    '--min-abundance', 
    type=int,
    default=1, 
    show_default=True,
    help='Minimum absolute abundance of a feature to plot.'
    )
def run_app(port, input, wait_time, relative, rm_unmapped, min_abundance):
   run_server(port, input, wait_time, relative, rm_unmapped, min_abundance) 

if __name__ == "__main__":
    cli()
