import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

# https://hackersandslackers.com/plotly-dash-with-flask/
def create_dashboard(server):
    """Create a Plotly Dash dashboard."""
    global dash_app

    dash_app = dash.Dash(
        server=server,
        routes_pathname_prefix='/dashapp/',
        external_stylesheets=[
            '/static/css/styles.css',
        ]
    )

    # Create Dash Layout
    dash_app.layout = html.Div(id='dash-container')
    init_callbacks(dash_app)

    return dash_app.server

def init_callbacks(dash_app):
    @dash_app.callback(Output('heatmap-id', 'figure'),
                       [Input('experiment-id', 'value')])
    def update_graph(vid):
        print("callback", vid)
        return vid
