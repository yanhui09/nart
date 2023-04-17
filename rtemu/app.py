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

@app.route('/')
def index():
    return render_template('index.html')

def run_server(port, input, wait_time):
    """
    Run the server.
    :param port: Port to run the app on.
    :param input: Path to the input TSV file.
    :param wait_time: Time to wait (in minutes) if input file is missing.
    """
    if not os.path.exists(input):
        print(f"Input file '{input}' not found. Waiting for {wait_time} minute(s)...")
        time.sleep(wait_time * 60)
        if not os.path.exists(input):
            print("I'm waiting to be fed. ;)")
            return

    @app.route('/plot')
    def plot():
        return generate_plot(input)

    app.run(debug=False, port=port)

    last_mod_time = os.stat(input).st_mtime
    while True:
        time.sleep(1)
        if os.path.exists(input):
            mod_time = os.stat(input).st_mtime
            if mod_time > last_mod_time:
                last_mod_time = mod_time
                with app.app_context():
                    generate_plot(input)
        else:
            print(f"Input file '{input}' not found. Waiting for {wait_time} minute(s)...")
            time.sleep(wait_time * 60)

if __name__ == '__main__':
    run_server()
