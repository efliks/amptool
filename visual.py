import json
import dash
import dash_cytoscape as cyto
import dash_html_components as html


app = dash.Dash(__name__)

atom_colors = {
    'H': 'cyan',
    'O': 'red',
    'N': 'blue',
    'C': 'gray'
}


def get_atom(node):
    return {
        'data': {
            'id': node['id'],
            'label': node['symbol'] + '-' + str(node['id'] + 1)
        },
        'classes': atom_colors[node['symbol']]
    }


def get_bond(node):
    return {
        'data': {
            'source': node['source'],
            'target': node['target']
        }
    }


stylesheet = [
    {
        'selector': 'node',
        'style': {
            'content': 'data(label)'
        }
    }
]
for _, color in atom_colors.items():
    stylesheet.append({
        'selector': '.' + color,
        'style': {
            'background-color': color
        }
    })


with open('molecule.json') as f:
    js = json.load(f)
    els = list(map(get_atom, js['nodes']))
    els.extend(list(map(get_bond, js['links'])))

app.layout = html.Div([
    cyto.Cytoscape(
        id='cytoscape-two-nodes',
        layout={'name': 'cose'},
        style={'width': '800px', 'height': '600px'},
        elements=els,
        stylesheet=stylesheet
    )
])

if __name__ == '__main__':
    app.run_server(debug=True)
