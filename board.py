import dash
import dash_auth
import dash_core_components as dcc
import dash_html_components as html
import dash_table as dt

#external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__)  # , external_stylesheets=external_stylesheets)
app.config.suppress_callback_exceptions = True

VALID_USERNAME_PASSWORD_PAIRS = {
    'hello': 'world'
}

auth = dash_auth.BasicAuth(
    app,
    VALID_USERNAME_PASSWORD_PAIRS
)

app.layout = html.Div(children=[
    html.H1(children='Hello Dash'),

    html.Div(children='''
        Dash: A web application framework for Python.
    '''),

    dcc.Graph(
        id='example-graph',
        figure={
            'data': [
                {'x': [1, 2, 3], 'y': [4, 1, 2], 'type': 'bar', 'name': 'SF'},
                {'x': [1, 2, 3], 'y': [2, 4, 5],
                    'type': 'bar', 'name': u'Montréal'},
            ],
            'layout': {
                'title': 'Dash Data Visualization'
            }
        }
    ),
    html.Div(children=html.Canvas(
        id="chemdoodle-viewer", width=386, height=386))
])

if __name__ == '__main__':
    app.run_server(debug=False, host='0.0.0.0', port=8080)

# <canvas id="chemdoodle-viewer" width="386" height="386" class="ChemDoodleWebComponent" style="width: 386px; height: 386px; background-color: rgb(255, 255, 255);"></canvas>
