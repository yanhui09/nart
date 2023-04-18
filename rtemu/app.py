import os
import time
import click
import pandas as pd
import plotly.graph_objs as go
from flask import Flask, render_template


app = Flask(__name__)
app.config['SECRET_KEY'] = 'secret!'

def generate_plot(input_file, relative=False, rm_unmapped=False):
    df = pd.read_csv(input_file, sep='\t', index_col=0)
    if rm_unmapped:
        df = df.drop('unassigned', axis=0, errors='ignore')
    if relative:
        df[df.columns[1:]] = df[df.columns[1:]].div(df[df.columns[1:]].sum(axis=0)) * 100
     
    sample_names = df.columns.tolist()[1:]
    # Get a list of all unique taxa in the dataframe
    taxa = df['taxonomy'].unique()

    data = []
    for i, taxon in enumerate(taxa):
        # Create a subset dataframe for the current taxon
        taxon_df = df[df['taxonomy'] == taxon]

        # Create a bar trace for each sample, with the abundance of the taxon in the sample
        trace = go.Bar(
            x=sample_names,
            y=taxon_df[sample_names].values[0],
            name=taxon,
            marker_color=f'rgba({(i*30)%255}, {(i*60)%255}, {(i*90)%255}, 0.7)'
        )
        data.append(trace)

    layout = go.Layout(
        barmode='stack',
        title='Taxonomy Composition',
        xaxis=dict(title='Barcode'),
        yaxis=dict(title='Abundance')
    )

    fig = go.Figure(data=data, layout=layout)
    return fig.to_json()

def get_modification_time(input_file):
    """
    Get the last modification time of the input file.
    """
    mod_time = os.stat(input_file).st_mtime
    mod_time_str = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(mod_time))
    return mod_time_str

def get_num_files(input_dir, file_type):
    """
    Get the number of files in the input directory.
    """
    return len([f for f in os.listdir(input_dir) if f.endswith(file_type)])

def get_latest_file(input_dir, file_type):
    """
    Get the latest file in the input directory.
    """
    files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith(file_type)]
    latest_file = max(files, key=os.path.getctime)
    return latest_file

def get_num_uniq_lines(input_file):
    """
    Get the number of unique lines in the input file.
    """
    with open(input_file) as f:
        lines = f.readlines()
    return len(set(lines))

@app.route('/')
def index():
    # use the file in /templates in rtemu package
    input_file = os.path.join(os.path.dirname(__file__), 'templates', 'input.tsv')
    #mod_time_str = get_modification_time(input_file)
    mod_time_str = '0000-00-00 00:00:00'
    plot_json = generate_plot(input_file, relative=False, rm_unmapped=False)
    num_fqs = 0
    num_batches = 0
    latest_fq = '0000-00-00 00:00:00'
    latest_batch = '0000-00-00 00:00:00'
    pct_complete = 0
    return render_template('index.html', 
                           plot_json=plot_json, mod_time_str=mod_time_str,
                           num_fqs=num_fqs, num_batches=num_batches, 
                           latest_fq=latest_fq, latest_batch=latest_batch,
                           pct_complete=pct_complete)

def run_server(port, work_dir, wait_time, relative, rm_unmapped):
    """
    Run the server.
    :param port: Port to run the app on.
    :param work_dir: Path to the work_dir.
    :param wait_time: Time to wait (in minutes) if input file is missing.
    :param relative: Whether to plot in relative abundance.
    :param rm_unmapped: Whether to remove unmapped reads.
    """
    input_file = os.path.join(work_dir, 'otu_table.tsv')
    fq_txt = os.path.join(work_dir, 'fqs.txt')
    batch_dir = os.path.join(work_dir, 'batches')

    if not os.path.exists(input_file):
        print(f"Input file '{input_file}' not found. Waiting for {wait_time} minute(s)...")
        time.sleep(wait_time * 60)
        if not os.path.exists(input_file):
            print("I'm waiting to be fed. ;)")
            return
        
    @app.route('/plot')
    def plot():
        return generate_plot(input_file, relative, rm_unmapped)
    
    @app.route('/mod_time')
    def mod_time():
        return get_modification_time(input_file)
    
    @app.route('/num_fqs')
    def num_fqs():
        return str(get_num_uniq_lines(fq_txt))
    
    @app.route('/latest_fq')
    def latest_fq():
        return get_modification_time(fq_txt)
    
    @app.route('/num_batches')
    def num_batches():
        return str(get_num_files(batch_dir, '.tsv'))
    
    @app.route('/latest_batch')
    def latest_batch():
        return str(get_modification_time(get_latest_file(batch_dir, '.tsv')))
    
    @app.route('/pct_complete')
    def pct_complete():
        return str(get_num_files(batch_dir, '.tsv') / get_num_uniq_lines(fq_txt) * 100)
    
    app.run(debug=False, port=port)

if __name__ == '__main__':
    run_server()
