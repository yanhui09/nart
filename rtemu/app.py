import os
import time
import click
import pandas as pd
import plotly.graph_objs as go
from flask import Flask, render_template


app = Flask(__name__)
app.config['SECRET_KEY'] = 'secret!'

def generate_plot(input_file):
    df = pd.read_csv(input_file, sep='\t', index_col=0)
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

@app.route('/')
def index():
    # use the file in /templates in rtemu package
    input_file = os.path.join(os.path.dirname(__file__), 'templates', 'input.tsv')
    mod_time_str = get_modification_time(input_file)
    plot_json = generate_plot(input_file)
    return render_template('index.html', plot_json=plot_json, mod_time_str=mod_time_str)

def run_server(port, input_file, wait_time):
    """
    Run the server.
    :param port: Port to run the app on.
    :param input_file: Path to the input TSV file.
    :param wait_time: Time to wait (in minutes) if input file is missing.
    """
    if not os.path.exists(input_file):
        print(f"Input file '{input_file}' not found. Waiting for {wait_time} minute(s)...")
        time.sleep(wait_time * 60)
        if not os.path.exists(input_file):
            print("I'm waiting to be fed. ;)")
            return

    @app.route('/plot')
    def plot():
        return generate_plot(input_file)
    
    @app.route('/mod_time')
    def mod_time():
        return get_modification_time(input_file)
    
    app.run(debug=False, port=port)

if __name__ == '__main__':
    run_server()
