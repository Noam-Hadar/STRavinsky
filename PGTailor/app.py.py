#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import dash
from dash.dependencies import Input, Output
from dash import dcc
from dash import html
from dash import dash_table
from itertools import product
import urllib.parse
import pandas as pd
import plotly.graph_objects as go
from flask import request


# In[ ]:


external_stylesheets = ['assets/PGTailor.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
server = app.server
app.title = "PGTailor"

chromosomes = [str(i) for i in range(1,23)] + ['X']
hg19_sizes = {
    '1' : 249250621,
    '2' : 243199373,
    '3' : 198022430,
    '4' : 191154276,
    '5' : 180915260,
    '6' : 171115067,
    '7' : 159138663,
    '8' : 146364022,
    '9' : 141213431,
    '10' : 135534747,
    '11' : 135006516,
    '12' : 133851895,
    '13' : 115169878,
    '14' : 107349540,
    '15' : 102531392,
    '16' : 90354753,
    '17' : 81195210,
    '18' : 78077248,
    '19' : 59128983,
    '20' : 63025520,
    '21' : 48129895,
    '22' : 51304566,
    'X' : 155270560,
    'Y' : 59373566
}

hg38_sizes = {
    '1' : 248956422,
    '2' : 242193529,
    '3' : 198295559,
    '4' : 191154276,
    '5' : 181538259,
    '6' : 170805979,
    '7' : 159345973,
    '8' : 145138636,
    '9' : 138394717,
    '10' : 133797422,
    '11' : 135086622,
    '12' : 133275309,
    '13' : 114364328,
    '14' : 107043718,
    '15' : 101991189,
    '16' : 90338345,
    '17' : 83257441,
    '18' : 80373285,
    '19' : 58617616,
    '20' : 64444167,
    '21' : 46709983,
    '22' : 50818468,
    'X' : 156040895,
    'Y' : 57227415
}

T2T_sizes = {
    '1' : 248387328,
    '2' : 242696752,
    '3' : 201105948,
    '4' : 193574945,
    '5' : 182045439,
    '6' : 172126628,
    '7' : 160567428,
    '8' : 146259331,
    '9' : 150617247,
    '10' : 134758134,
    '11' : 135127769,
    '12' : 133324548,
    '13' : 113566686,
    '14' : 101161492,
    '15' : 99753195,
    '16' : 96330374,
    '17' : 84276897,
    '18' : 80542538,
    '19' : 61707364,
    '20' : 66210255,
    '21' : 45090682,
    '22' : 51324926,
    'X' : 154259566,
    'Y' : 62460029
}


chromosomalSizes = {'hg19' : hg19_sizes, 'hg38' : hg38_sizes, 'T2T' : T2T_sizes}


# In[ ]:


app.layout = html.Div([
    html.Div([
        html.Img(
            src = "assets/BirkLab_logo.png",
            style={'width':'18%', 'float' : 'right'}
        ),
        html.Img(
            src = "assets/PGTailor.png",
            style={'width':'18%', 'float' : 'left'}
        ),
        ],
        style={'width': '100%', 'display' : 'inline-block','marginLeft' : 'auto', 'marginRight' : 'auto' , 'textAlign' : 'center'}
    ),
    html.Div(
    [
        html.Div(html.H1('Parameters'), style = {'font-family' : 'gisha', 'width': '100%', 'display': 'inline-block','marginLeft':'auto', 'marginRight':'auto', 'textAlign' : 'center'}),
        html.Div([
            dcc.RadioItems(
            id = 'referenceRadioButtons',
            options=[
                {'label': 'hg38', 'value': 'hg38'},
                {'label': 'hg19', 'value': 'hg19'},
                {'label': 'T2T', 'value': 'T2T'},
            ],
            value='hg38',
            labelStyle={'display': 'inline-block'},
            style={'display': 'inline-block','marginLeft':'auto', 'marginRight':'auto', 'textAlign' : 'left'}
        ),
        html.Span(html.A('Convert coordinates', href = 'https://liftover.broadinstitute.org/', target='_blank', style = {'display' : 'inline-block'}), style = {'display' : 'inline-block', 'marginLeft' : 20}),
        html.Div([
            html.Div(
            dcc.Dropdown(
                id = 'chromosome',
                options = [{'label': 'chr' + str(v), 'value': str(v)} for v in chromosomes],
                multi = False,
                placeholder = 'Chromosome',
                style = {'textAlign' : 'center', 'font-family' : 'gisha'}
            ), style = {'width' : '28%', 'display' : 'inline-block', 'transform' : 'translateY(35%)'}),
            html.Div(
                html.P(':', style = {'textAlign' : 'center', 'font-family' : 'gisha'}),
                style = {'width' : '5%', 'display' : 'inline-block'}
            ),
            html.Div(
            dcc.Input(
                id = 'position',
                style = {'textAlign' : 'center', 'font-family' : 'gisha'},
                placeholder = 'Position',
                type = 'number',
                min = 1,
                max = 250000000,
                debounce = True,
                step = 1,
            ), style = {'width' : '60%','display' : 'inline-block'}),
            ],
            style = {'width' : '100%', 'display' : 'inline-block'}
        ),
        ]
        ),
        html.P(),
        html.Div([
            html.Span('Look at '),
            dcc.Input(
                id = 'flank_size',
                style = {'display' : 'inline-block','textAlign' : 'center'},
                type = 'number',
                min = 1,
                max = 3,
                value = 2,
                step = 1,
            ),
            html.Span('mbp from each side of the given coordinate.'),
        ]),
        html.Div([
            html.P(),
            html.Div([
                html.Span('4bp STRs must repeat at least '),
                dcc.Input(
                    id = '4bprepeats',
                    style = {'display' : 'inline-block','textAlign' : 'center'},
                    type = 'number',
                    min = 5,
                    max = 30,
                    value = 6,
                ),
                html.Span(' times.'),
            ]),
            html.P(),
            html.Div([
                html.Span('3bp STRs must repeat at least '),
                dcc.Input(
                    id = '3bprepeats',
                    style = {'display' : 'inline-block','textAlign' : 'center'},
                    type = 'number',
                    min = 5,
                    max = 30,
                    value = 11,
                ),
                html.Span(' times.'),
            ]),
            html.P(),
            html.Div([
                html.Span('2bp STRs must repeat at least '),
                dcc.Input(
                    id = '2bprepeats',
                    style = {'display' : 'inline-block','textAlign' : 'center'},
                    type = 'number',
                    min = 5,
                    max = 30,
                    value = 15,
                ),
                html.Span(' times.'),
            ]),
        ], style = {'display' : 'none'}),
        html.Div([
            html.P(),
            html.P(),
            html.Button(
                'Find STRs',
                id = 'first_analysis_button',
                disabled = True,
                n_clicks = 0,
                style = {'text-transform' : 'none', 'font-size' : '30px', 'width' : '100%'},
            ),
            ],
        style = {'width' : '100%', 'display' : 'inline-block','textAlign' : 'left'}
        )
    ],
    id = 'parameters_div',
    style={'width': '50%', 'display' : 'inline-block', 'textAlign' : 'left', 'float' : 'center'}
    ),
    
    html.Div(id  = 'str_select_div'),
    html.Div(id = 'chosen_strs_div'),
    
    html.Div(
        dcc.Loading(
            id = "loading_wheel",
            children = [html.Div([html.Div(id = "loading-output")])],
            type = "circle",
        ),
        style = {'width' : '100%', 'display' : 'inline-block','textAlign' : 'center', 'align' : 'center'},
        id = "loading_div"
    ),
    
    dcc.Store('reference_store'),
    dcc.Store('chromosome_store'),
    dcc.Store('position_store'),
    dcc.Store('flank_size_store'),
    dcc.Store('fourbprepeats'),
    dcc.Store('threebprepeats'),
    dcc.Store('twobprepeats'),
    dcc.Store('selected_STRs'),
    dcc.Store('genes'),
    dcc.Store('sequence'),
],
style={'width': '100%', 'display' : 'inline-block','textAlign' : 'center'}
   
)


# In[ ]:


@app.callback(Output('referenceRadioButtons', 'labelstyle'), [Input('2bprepeats', 'style')])
def get_ip(value):
    print(request.remote_addr)
    return {'display': 'inline-block'}


# In[ ]:


@app.callback(
    [Output('first_analysis_button', 'disabled'),
     Output('reference_store', 'data'),
     Output('chromosome_store', 'data'),
     Output('position_store', 'data'),
     Output('flank_size_store', 'data'),
     Output('fourbprepeats', 'data'),
     Output('threebprepeats', 'data'),
     Output('twobprepeats', 'data')
    ],
    [Input('referenceRadioButtons', 'value'),
     Input('chromosome', 'value'),
     Input('position', 'value'),
     Input('flank_size', 'value'),
     Input('4bprepeats', 'value'),
     Input('3bprepeats', 'value'),
     Input('2bprepeats', 'value')
    ]
)
def searchButtonAvailabilityStatus(reference_genome, chromosome, position, flank_size, fourbprepeats, threebprepeats, twobprepeats):
    disabled_btn = [True, reference_genome, chromosome, position, flank_size, fourbprepeats, threebprepeats, twobprepeats]
    enabled_btn = [False, reference_genome, chromosome, position, flank_size, fourbprepeats, threebprepeats, twobprepeats]
    if None in (reference_genome, chromosome, position, flank_size, fourbprepeats, threebprepeats, twobprepeats):
        return disabled_btn
    elif chromosomalSizes[reference_genome][chromosome] < int(position):
        return disabled_btn
    else:
        return enabled_btn


# In[ ]:


def designToolTip(r):
    toolTip = """
| Primer       | Sequence     | Tm     |
| :------------- | :----------: | -----------: |
| Forward   | fps    | fptm |
| Reverse   | rps | rptm |
| Nested   | nps | nptm |
"""
    if r['Forward primer'] == '':
        return 'Primer design did not succeed'
    toolTip = toolTip.replace('fps', r['Forward primer'])
    toolTip = toolTip.replace('fptm', r['Forward primer Tm'] + '°c')
    toolTip = toolTip.replace('rps', r['Reverse primer'])
    toolTip = toolTip.replace('rptm', r['Reverse primer Tm'] + '°c')
    toolTip = toolTip.replace('nps', r['Nested PCR primer'])
    toolTip = toolTip.replace('nptm', r['Nested PCR primer Tm'] + '°c')
    if 'failed' in toolTip:
        return 'Primer design did not succeed'
    if r['Sequence'] == 'Your mutation':
        return 'Your mutation'
    return toolTip

@app.callback(
    [Output('parameters_div', 'style'),
     Output('str_select_div', 'children'),
     Output('loading_wheel', 'children'),
     Output('genes', 'data')
    ],
    [Input('first_analysis_button', 'n_clicks'),
     Input('reference_store', 'data'),
     Input('chromosome_store', 'data'),
     Input('position_store', 'data'),
     Input('flank_size_store', 'data'),
     Input('fourbprepeats', 'data'),
     Input('threebprepeats', 'data'),
     Input('twobprepeats', 'data')
    ]
)
def startFirstAnalysis(n_clicks, reference, chromosome, position, flank_size, minFour, minThree, minTwo):
    print(n_clicks) #shit
    if n_clicks in [0, None]:
        return [{'width': '50%', 'display' : 'inline-block', 'textAlign' : 'left', 'float' : 'center'}, [], None, None]
    else:
        genes_df = pd.read_csv('assets/' + reference.lower() + '_map.csv')
        genes_df['Chr'] = genes_df['Chr'].astype(str)
        genes = genes_df[(genes_df['Chr'] == chromosome) & (genes_df['Start'] <= position + (flank_size * 1000000)) & (genes_df['End'] >= position - (flank_size * 1000000))].to_dict('records')
        df = pd.read_csv('assets/' + reference.replace('T2T', 'CHM13') + '_chr' + chromosome + '_STRs.csv.gz', dtype = 'str')
        df['start'] = df['Coordinates'].apply(lambda x : int(x.split(':')[1].split('-')[0]))
        df = df[(df['start'] > position - (flank_size * 1000000)) & (df['start'] < position + (flank_size * 1000000))]
        df['STR_len'] = df['Sequence'].str.split(' x ').str[0].str.len()
        df['repeat#'] = df['Sequence'].apply(lambda x : int(x.split(' x ')[1]))
        def checkIfEnoughRepeats(STR_len, repeatNumber):
            if STR_len == 4:
                return repeatNumber >= minFour
            elif STR_len == 3:
                return repeatNumber >= minThree
            else:
                return repeatNumber >= minTwo
        df.fillna('', inplace = True)
        df = df[df.apply(lambda x : checkIfEnoughRepeats(x['STR_len'], x['repeat#']),axis = 1)]    
        del df['STR_len']
        del df['repeat#']
        df['Distance'] = df['start'] - position
        del df['start']
        df['Primer lengths'] = df.apply(lambda x : ', '.join([str(len(x[i])) for i in ['Forward primer', 'Reverse primer', 'Nested PCR primer']]).replace('0, 0, 0','') , axis = 1)
        df['PrimerBLAST'] = df.apply(lambda x : '❌' if x['Primer lengths'] == '' else '[Check specificity🧬](https://www.ncbi.nlm.nih.gov/tools/primer-blast/?PRIMER_LEFT_INPUT=' + x['Forward primer'] + '&PRIMER_RIGHT_INPUT=' + x['Reverse primer'] + '&PRIMER_SPECIFICITY_DATABASE=PRIMERDB/genome_selected_species)', axis = 1)
        cols = ['Distance', 'PrimerBLAST', 'Sequence', 'Coordinates']
        others = [col for col in df.columns if col not in cols]
        df = df[cols + others]
        
        ddt = dash_table.DataTable(
            id = 'table',
            data = df.to_dict('records'),
            columns = [{'id': c, 'name': c, 'presentation': 'markdown'} if c == 'PrimerBLAST' else {'id': c, 'name': c} for c in df.columns if c in cols],
            page_action = 'none',
            row_selectable='multi',
            fixed_rows={'headers' : True},
            style_cell = {
                'textAlign': 'left',
                'overflow': 'hidden',
                'textOverflow': 'ellipsis',
                'maxWidth': 0
            },
            style_data_conditional=[
                {
                'if': {'row_index': 'odd'},
                'backgroundColor': 'rgb(247, 247, 255)'
                },
                {
                'if': {'row_index': 'even'},
                'backgroundColor': 'rgb(252, 252, 255)'
                },
                {
                'if': {'filter_query': '{Distance} = 0'},
                'backgroundColor': 'rgb(255, 200, 200)'
                },
                {
                "if": {"state": "selected"},
                "backgroundColor": "inherit !important",
                "border": "inherit !important",
                },
            ],
            style_data={
                'whiteSpace': 'normal',
                'height': 'auto',
                },
            style_header = {
                'backgroundColor': 'rgb(220, 220, 255)',
                'fontWeight': 'bold',
                'whiteSpace' : 'normal'
                },
            style_table={
                'maxHeight': '200px',
                'maxwidth' : '150%',
                'border': 'thin lightgrey solid'
            },
            #tooltip_data=[
            #    {
            #        column: {'value': str(value), 'type': 'markdown'}
            #        for column, value in row.items()
            #    } for row in df.to_dict('rows')
            #],
            tooltip_data=[
                {
                    column: {'value': designToolTip(row), 'type': 'markdown'}
                    for column, value in row.items()
                } for row in df.to_dict('rows')
            ],
            css=[{
                'selector': '.dash-table-tooltip',
                'rule': 'background-color: white; font-family: monospace; color: black'
            }],
            style_as_list_view = False
        )
        return [{'display' : 'none'}, [ddt], None, genes]


# In[ ]:


@app.callback(
    [Output('chosen_strs_div', 'children'),
     Output('selected_STRs', 'data')
    ],
    [Input('table', 'derived_virtual_selected_rows'),
     Input('table', 'derived_virtual_data'),
     Input('genes', 'data'),
     Input('position_store', 'data')
    ]
)
def manageChosen(selected_rows, data, genes, position):
    chosen_strs_div = []
    if selected_rows == None:
        return [None, None]
    selected_data = [data[row] for row in selected_rows]
    distances = [row['Distance'] for row in selected_data]
    before = len([n for n in distances if n < 0])
    after = len([n for n in distances if n > 0])
    before_and_after_message = 'You chose ' + str(before) + ' before and ' + str(after) + ' after the selected location'
    chosen_strs_div.append(html.P(before_and_after_message))
    data_df = pd.DataFrame(data)
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=[data_df['Distance'].min(), data_df['Distance'].max()], y=[0,0], mode='markers', text = ['left border', 'right border'], marker_size=10, marker=dict(color='black')))
    fig.add_trace(go.Scatter(x=[0], y=[0], mode='markers', text = ['Chosen position'], marker_symbol = 'star', marker_size = 20, marker=dict(color='red', opacity = 0.6)))
    if len(selected_data) > 0:
        df = pd.DataFrame(selected_data)
        X = df['Distance'].tolist()
        Y = [0 for y in X]
        fig.add_trace(go.Scatter(x=X, y=Y, mode='markers', text = df['Coordinates'] + ', ' + df['Sequence'], marker_size = 15, marker_symbol = 'diamond-tall', marker=dict(color='blue', opacity = 0.5)))
    for n, gene in enumerate(genes):
        start = gene['Start'] - position
        end = gene['End'] - position
        m = n%12
        fig.add_trace(go.Scatter(x = [start, end], y = [-5 * (m + 1), -5 * (m + 1)], mode = 'lines', text = gene['Gene']))
        fig.add_trace(go.Scatter(x = [(start + end)/2], y = [-5 * (m + 1) + 1], mode = 'text', text = gene['Gene']))
    fig.update_xaxes(showgrid = False)
    fig.update_yaxes(showgrid = False, zeroline = True, zerolinecolor='black', zerolinewidth=3, showticklabels=False)
    fig.update_layout(height = 600, plot_bgcolor='white')
    fig.update_layout(showlegend=False) 
    chosen_strs_div.append(dcc.Graph(figure = go.Figure(fig)))
    if len(selected_data) > 0:
        col_order = ['Distance', 'Sequence', 'Coordinates', 'Forward primer', 'Forward primer Tm', 'Reverse primer', 'Reverse primer Tm', 'Product size', 'Nested PCR primer', 'Nested PCR primer Tm', 'Nested PCR primer orientation', 'Nested PCR product size']
        csv_report = pd.DataFrame(selected_data)[col_order]
        csv_report = csv_report.fillna('').to_csv(na_rep='', index = False).replace(',nan,', ',,')
        download_href = "data:text/csv;charset=utf-8," + urllib.parse.quote(csv_report)
        linkStyle = {'marginRight' : 20, 'display' : 'inline-block', 'color' : "#077be2",'font-family' : 'gisha', 'fontSize' : 15, 'textAlign' : 'center' ,'borderWidth' : '2px','borderColor' : "#000044",'borderStyle' : 'groove','padding': 4,'borderRadius' : '5px'}
        chosen_strs_div = [chosen_strs_div[0], html.A('Download chosen STRs', href = download_href, download = 'chosen STRs.csv', style = linkStyle), chosen_strs_div[-1]]
    else:
        chosen_strs_div = [chosen_strs_div[0], html.Br(),html.Br(), chosen_strs_div[-1]]
    return [chosen_strs_div, selected_data]


# In[ ]:


if __name__ == '__main__':
    #app.run_server(port = 2156)
    server.run(debug=True)

