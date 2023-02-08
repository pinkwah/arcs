import dash
#import dash_core_components as dcc
#import dash_html_components as html
import dash_bootstrap_components as dbc
from dash_bootstrap_templates import load_figure_template
from dash import html
from dash import dash_table
from dash import dcc
from dash import ctx
from dash.dependencies import Input, Output, State
import plotly.express as px
import pandas as pd
import time
import base64
import io
import numpy as np 
from monty.serialization import loadfn
from arcs.setup_functions import GenerateInitialConcentrations
from arcs.analysis_2 import AnalyseSampling
from arcs.traversal import Traversal
import pickle
import warnings


external_stylesheets = [
     dbc.themes.QUARTZ, 
]

load_figure_template('QUARTZ')

######data and sliders

def load_data(filename): # has to be a .json file in dict format
    data = loadfn(filename)
    return(data)

def keys_by_depth(dict_, depth=0, output=None):
    if output is None:
        output = {}
    if not depth in output:
        output[depth] = set()
    for key in dict_:
        output[depth].add(key)
        if isinstance(dict_[key], dict):
            keys_by_depth(dict_[key], depth+1, output)
    return output  

def make_sliders(import_data,labels):
    
    def sliderform(key,slider,label):
        # this is very janky - need to find a better way of storing the data
        correct_order = sorted([float(x) for x in list(slider)])
        #marks = {float(x):str('{:.2F}'.format(float(x))) for x in list(correct_order)}
        marks = {int(x):str(x) for x in correct_order}
        minval = float(list(correct_order)[0])
        maxval = float(list(correct_order)[-1])
        
        return(html.Div(style={'display':'flex',
                              'padding':'0.8rem'},
                        children=[
                            dbc.Col(
                                children=html.Label(children=label[key]),
                                width=2
                            ),
                            dbc.Col(
                                children=dcc.Slider(id='slider-{}'.format(key),
                                                    min=minval,
                                                    max=maxval,
                                                    step=None,
                                                    marks=marks,
                                                    value=minval,
                                                    updatemode='drag',
                                                    tooltip={"placement": "bottom", "always_visible":True}),
                                width=10,                         
                            )
                        ]
                       )
              )  
    
    slider_keys = keys_by_depth(import_data)  
    sliders = []
    for k in slider_keys:
        sliders.append(sliderform(k,sorted(slider_keys[k]),labels))
    
    return(sliders)


#run data fields
g = pickle.load(open('./assets/data/SCAN_graph.p','rb'))
t = Traversal(graph=g,reactions='./assets/data/SCAN_reactions.p')

backgroundcolours='rgba(0,0,0,0.1)'
graph = None
table = None
table2 = None
table3 = None
table4 = None
table5 = None
meta = None
sliders = make_sliders(g,labels={0:'T (K)',1:'P (bar)'})
gic = GenerateInitialConcentrations(g)
gic.all_zero(include_co2=False)
gic.update_ic({'SO2':10e-6,'NO2':50e-6,'H2S':30e-6,'H2O':20e-6})
concs = gic.ic
settings = {'nprocs':1,
       'sample_length':320,
       'max_rank':10,
       'max_compounds':5,
       'probability_threshold':0.1,
       'path_depth':5,
       'ceiling':2000,
       'scale_highest':0.2}
ambient_settings = {'T':None,'P':None}










######################starting#######################
app = dash.Dash(__name__,external_stylesheets=external_stylesheets)

app.layout = html.Div(
    style={'padding':'2rem'},
    children=[
    dcc.Loading(
        id="loading-1", 
        type='cube',
        color='lightblue',
        fullscreen=True,
        children=html.Div(id="loading-output-1"),
        style={'background-color':'rgba(0.1,0.1,0.1,0.5)'}
        ),
    dbc.Row(
        className='Row',
        style={'margin':'1rem','display':'flex','justify-content': 'space-between','flex-wrap': 'wrap'},
        children=[
            #1st column
            dbc.Col(style={'align':'flex',
                            'boxShadow': '2px 2px 2px 2px pink', 
                            'border-radius': '5px',
                            'padding':'1rem',
                            'backgroundColor':backgroundcolours},
                     children=[
                         html.H5('Control'),
                         dbc.Row(
                             children=sliders
                         ),
                         dbc.Row(
                             style={'padding':'2rem',
                                   'align':'center'},
                             children = [
                                 dbc.Col(
                                     children = [
                                         html.H6('Concentrations (PPM)'),
                                         html.Div(
                                             style={'padding':'0.5rem'},
                                             children=dcc.Textarea(id='concs-out',
                                                      value="SO2=10,\nNO2=50,\nH2S=30,\nH2O=20",
                                                      style={
                                                          'backgroundColor':backgroundcolours,
                                                          'width':'10rem',
                                                          'height':'21rem',
                                                          'color':'white',
                                                      'border':'2px white solid'})),
                                         dbc.Button(style={'border':'2px white solid',
                                                      'marginTop':20,
                                                      'marginLeft':10},
                                                    children='update',
                                                    id='concs-submit',
                                                    n_clicks=0,
                                                    className="me-1",
                                                    color='Primary')

                                     ]
                                 ),
                                 #########################settings
                                 dbc.Col(
                                     style={'align':'center'},
                                     children = [
                                         html.H6('Algorithm Settings'),
                                         html.Div(
                                             style={'padding':'0.5rem'},
                                             children=[html.P('samples:'),
                                                       dcc.Input(id='samples',
                                                                 value="10",
                                                                 debounce=True,
                                                                 style={'backgroundColor':backgroundcolours,
                                                                        'color':'white',
                                                                        'border':'2px white solid'}
                                                                )]),
                                         html.Div(
                                             style={'padding':'0.5rem'},
                                             children=[html.P('path depth:'),
                                                       dcc.Input(id='pathdepth',
                                                                 value="5",
                                                                 debounce=True,
                                                                 style={'backgroundColor':backgroundcolours,
                                                                        'color':'white',
                                                                        'border':'2px white solid'}
                                                                )]),
                                         html.Div(
                                             style={'padding':'0.5rem'},
                                             children=[html.P('probability cutoff (max=1):'),
                                                       dcc.Input(id='probability_cutoff',
                                                                 value="0.01",
                                                                 debounce=True,
                                                                 style={'backgroundColor':backgroundcolours,
                                                                        'color':'white',
                                                                        'border':'2px white solid'}
  
                                                                )]), 
                                         html.Div(
                                             style={'padding':'0.5rem'},
                                             children=[html.P('Ceiling (%):'),
                                                       dcc.Input(id='ceiling',
                                                                 value="2000",
                                                                 debounce=True,
                                                                 style={'backgroundColor':backgroundcolours,
                                                                        'color':'white',
                                                                        'border':'2px white solid'}
                                                                )]), 
                                         
                                         html.Div(
                                             style={'padding':'0.5rem'},
                                             children=[html.P('scale large concentrations:'),
                                                       dcc.Input(id='scale_highest',
                                                                 value="0.1",
                                                                 debounce=True,
                                                                 style={'backgroundColor':backgroundcolours,
                                                                        'color':'white',
                                                                        'border':'2px white solid'}
                                                                )]), 
                                     ]
                                 ),
                                 dbc.Col(
                                     style={'align':'center'},
                                     children = [
                                         html.H6(' '),
                                         html.Div(
                                             style={'padding':'0.5rem'},
                                             children=[html.P('max rank:'),
                                                       dcc.Input(id='max_rank',
                                                                 value="10",
                                                                 debounce=True,
                                                                 style={'backgroundColor':backgroundcolours,
                                                                        'color':'white',
                                                                        'border':'2px white solid'}
                                                                )]),
                                         html.Div(
                                             style={'padding':'0.5rem'},
                                             children=[html.P('max compounds:'),
                                                       dcc.Input(id='max_compounds',
                                                                 value="5",
                                                                 debounce=True,
                                                                 style={'backgroundColor':backgroundcolours,
                                                                        'color':'white',
                                                                        'border':'2px white solid'}
                                                                )]),    
                                         html.Div(
                                             style={'padding':'0.5rem','width':'15rem'},
                                             children=[html.P('shortest path algorithm:'),
                                                       dcc.Dropdown(
                                                           id='method',
                                                          # style={'backgroundColor':backgroundcolours,
                                                          #              'color':'white',
                                                          #              'border':'2px white solid'},
                                                           options={'Bellman-Ford': 'Bellman-Ford',
                                                               'Dijkstra': 'Dijkstra'},
                                                           value='Bellman-Ford'),
                                                      ]),  
                                         html.Div(
                                             style={'padding':'0.5rem','width':'15rem'},
                                             children=[html.P(['include CO',html.Sub(2),':']),
                                                       dcc.Dropdown(
                                                           id='include_co2',
                                                           #style={'backgroundColor':backgroundcolours,
                                                           #             'color':'white',
                                                           #      },#'border':'2px white solid'},
                                                           options={'True': 'True',
                                                               'False': 'False'},
                                                           value='False',
                                                           className="me-1")
                                                                ]),  
                                         html.Div(
                                             children=dbc.Button(style={'border':'2px white solid',
                                                                        'marginTop':30},
                                                                 children='Submit',
                                                                 id='submit-val',
                                                                 n_clicks=0,
                                                                 className="me-1",
                                                                 color='Primary')
                                         )
                                     ]
                                 )
                             ]
                         )
                     ],
                    
                   ),
            #1st column
            dbc.Col(
                    style={'align':'center',
                           'boxShadow': '2px 2px 2px 2px pink', 
                           'border-radius': '5px',
                           'padding':'1rem', 
                           'backgroundColor':backgroundcolours},
                    children=[
                        html.H5('Change in Concentrations'),
                        html.Div(id='output-graph',
                                 children=graph
                                ),
                        html.Div(id='placeholder1', 
                                 children=None,
                                 style={'display':'none'}),
                        html.Div(id='placeholder2', 
                                 children=None,
                                 style={'display':'none'}),
                        html.Div(id='placeholder3', 
                                 children=None,
                                 style={'display':'none'}),
                        dcc.Loading(id='loading1',
                                    type='default',
                                    children=html.Div(id='loading-output'))
                             ],
                width={'offset':1}
                
                   ),
        ],
    ),
        dbc.Row(
            style={'margin':'1rem','display':'flex','justify-content': 'space-between','flex-wrap': 'wrap'},
            children=[
            dbc.Col(style={'align':'center',
                           'boxShadow': '2px 2px 2px 2px pink',  
                           'border-radius': '5px',
                           'padding':'1rem',
                           'backgroundColor':backgroundcolours},
                    children=[
                        html.H5('Most Frequent Reactions'),
                        html.Div(
                            id = 'reaction-stats',
                            style={'align':'center'},
                            children=table,
                            className='table'
                        )
                    ]
                   ),
            dbc.Col(style={'align':'center',
                           'boxShadow': '2px 2px 2px 2px pink', 
                           'border-radius': '5px',
                           'padding':'1rem',
                           'backgroundColor':backgroundcolours},
                    children=[
                        html.H5('Most Frequent Paths'),
                        html.Div(id = 'reaction-paths',
                                 style={'align':'center'},
                                 children=table2,
                                 className='table'
                        )
                    ],
                    width={'offset':1})
        ],   
    ),
        dbc.Row(
            style={'margin':'1rem','display':'flex','justify-content': 'space-between','flex-wrap': 'wrap'},
            children=[
                dbc.Col(style={'align':'center',
                               'boxShadow': '2px 2px 2px 2px pink',  
                               'border-radius': '5px',
                               'padding':'1rem',
                               'backgroundColor':backgroundcolours},
                    children=[
                        html.H5('Meta Data'),
                        html.Div(
                            id = 'metadata',
                            style={'align':'center'},
                            children=table3,
                            #className='table'
                        )
                    ],
                   ),
                dbc.Col(style={'align':'center',
                               'boxShadow': '2px 2px 2px 2px pink',  
                               'border-radius': '5px',
                               'padding':'1rem',
                               'backgroundColor':backgroundcolours},
                    children=[
                        html.H5('Initial Concentrations:'),
                        html.Div(
                            id = 'init_concs',
                            style={'align':'center'},
                            children=table4,
                            #className='table'
                        )
                    ],
                        width={'offset':1}
                   ),
                dbc.Col(style={'align':'center',
                               'boxShadow': '2px 2px 2px 2px pink',  
                               'border-radius': '5px',
                               'padding':'1rem',
                               'backgroundColor':backgroundcolours},
                    children=[
                        html.H5('Final Concentrations:'),
                        html.Div(
                            id = 'final_concs',
                            style={'align':'center'},
                            children=table5,
                            #className='table'
                        )
                    ],
                        width={'offset':1}
                   ),
            ]
        ),               

        dbc.Row(
            children=[
                dbc.Col(
                    children=[
                        html.Div(style={'display':'inline-block','align':'end'},
                                 children=[
                                     html.Img(src=app.get_asset_url('images/logos.png'),
                                              style={'width':'50%','height':'50%','padding':'0.5rem','align':'end'}),
                                 ]
                                ),
                    ]
                )

            ]
        )
    
])
######### update concentrations 
@app.callback(Output('placeholder1','children'),
              Input('concs-submit','n_clicks'),
              State('concs-out','value'))    

def update_concentrations(n_clicks,input1):
    if n_clicks > 0:
        itext = str(input1).split(',')
        for i in itext:
            try:
                spec,num = i.split('=')
                if '\n' in spec:
                    spec = spec.split('\n')[1]
                concs[spec] = float(num)*1e-6
            except:
                pass
            
############ update settings            
@app.callback(Output('placeholder2','children'),
              [Input('samples','value'),
               Input('pathdepth','value'),
               Input('probability_cutoff','value'),
               Input('ceiling','value'),
               Input('scale_highest','value'),
               Input('max_rank','value'),
               Input('max_compounds','value'),
               Input('method','value'),
               Input('include_co2','value')])

def update_settings(*inputs):
    settings['sample_length'] = int(inputs[0])
    settings['path_depth'] = int(inputs[1])
    settings['probability_threshold'] = float(inputs[2])
    settings['ceiling'] = int(inputs[3])
    settings['scale_highest'] = float(inputs[4])
    settings['max_rank'] = int(inputs[5])
    settings['max_compounds'] = int(inputs[6])
    settings['method'] = inputs[7]
    settings['include_co2'] = bool(inputs[8])
    
########## update T and P     
@app.callback(Output('placeholder3','children'),
              [Input('slider-0','value'),
               Input('slider-1','value')]) 
              
def update_t_and_p(*inputs):
    ambient_settings['T'] = int(inputs[0])
    ambient_settings['P'] = int(inputs[1])
    
    
@app.callback([Output('metadata','children'),
               Output('reaction-stats','children'),
               Output('reaction-paths','children'),
               Output('init_concs','children'),
               Output('final_concs','children'),
               Output('output-graph','children'),
               Output('loading-output-1','children')
              ],
              Input('submit-val','n_clicks'))  

def apprun(btn1):
    if 'submit-val' == ctx.triggered_id:
        warnings.simplefilter('ignore')
        t.run(trange=[ambient_settings['T']],
              prange=[ambient_settings['P']],
              save=False,
              ic=concs,
              **settings)
        
        metadata = t.metadata
        #del metadata['init_concs']
        metadata = pd.Series(metadata).reset_index()
        metadata = metadata.rename(columns={0:'value'})
        
        metadata_table = dash_table.DataTable(
            columns=[{'name':i,'id':i} for i in metadata.columns],
            data=metadata.to_dict('records'),
            style_as_list_view=True,
            cell_selectable=False,
            style_cell={'font_family':'helvetica',
                    'text_align':'center',
                   'hover-background-color': '#555555'},
            style_header={'backgroundColor': 'Transparent'},
            style_data={'backgroundColor': 'Transparent',
                    'border':'none'},
            style_table={'height': '400px', 
                     'overflowY': 'auto'}
                                    )
        
        
        init_concs = pd.Series(t.concs)[pd.Series(t.concs) > 0].to_dict()
        ic = {k:float('{:.2F}'.format(v/1e-6)) for k,v in init_concs.items()}
        
        del ic['CO2']
        df_i = pd.Series(ic)
        df_i = df_i.reset_index()
        df_i = df_i.T.rename({'index':'compound',0:'concentration (ppm)'})
        df_i = df_i.T
        
        init_concs_table = dash_table.DataTable(
            columns=[{'name':i,'id':i} for i in df_i.columns],
            data=df_i.to_dict('records'),
            style_as_list_view=True,
            cell_selectable=False,
            style_cell={'font_family':'helvetica',
                    'text_align':'center',
                   'hover-background-color': '#555555'},
            style_header={'backgroundColor': 'Transparent'},
            style_data={'backgroundColor': 'Transparent',
                    'border':'none'},
            style_table={'height': '400px', 
                     'overflowY': 'auto'}
        )        

        fc_i = pd.Series(t.final_concs[ambient_settings['T']][ambient_settings['P']])[pd.Series(t.final_concs[ambient_settings['T']][ambient_settings['P']]) >0].to_dict()
        fc = {k:float('{:.2F}'.format(v)) for k,v in fc_i.items()}
        df_f = pd.Series(fc)
        df_f = df_f.reset_index()
        df_f = df_f.T.rename({'index':'compound',0:'concentration (ppm)'})
        df_f = df_f.T
        
        final_concs_table = dash_table.DataTable(
            columns=[{'name':i,'id':i} for i in df_f.columns],
            data=df_f.to_dict('records'),
            style_as_list_view=True,
            cell_selectable=False,
            style_cell={'font_family':'helvetica',
                    'text_align':'center',
                   'hover-background-color': '#555555'},
            style_header={'backgroundColor': 'Transparent'},
            style_data={'backgroundColor': 'Transparent',
                    'border':'none'},
            style_table={'height': '400px', 
                     'overflowY': 'auto'}
                                    )
        
        analyse =  AnalyseSampling(t.data)
        analyse.mean_sampling() ; mean = analyse.mean_data[ambient_settings['T']][ambient_settings['P']]
        analyse.reaction_statistics() ; stats = pd.DataFrame(analyse.stats[ambient_settings['T']][ambient_settings['P']])
        analyse.reaction_paths() ; paths = pd.DataFrame(analyse.common_paths[ambient_settings['T']][ambient_settings['P']])


        stats_table = dash_table.DataTable(
            columns=[{'name':i,'id':i} for i in stats.columns],
            data=stats.to_dict('records'),
            style_as_list_view=True,
            cell_selectable=False,
            style_cell={'font_family':'helvetica',
                        'text_align':'center',
                       'hover-background-color': '#555555'},
            style_header={'backgroundColor': 'Transparent'},
            style_data={'backgroundColor': 'Transparent',
                        'border':'none'},
            style_table={'height': '400px',
                        'overflowY':'auto'},
        )
        
        
        paths_table = dash_table.DataTable(
            columns = [{'name':'paths','id':'paths'},{'name':'k','id':'k'},{'name':'frequency','id':'frequency'}],
        #columns=[{'name':i,'id':i} for i in new_stats.columns],
            data=paths.to_dict('records'),
            style_as_list_view=True,
            cell_selectable=False,
            style_cell={'font_family':'helvetica',
                        'text_align':'center',
                        'whiteSpace': 'pre-line',
                        'padding': '5px'},
            style_header={'backgroundColor': 'Transparent'},
            style_data={'backgroundColor': 'Transparent',
                        'border':'none'},
            style_table={'height': '400px', 
                         'overflowY': 'auto'}
                                    )
        df_m = pd.DataFrame({'comps':list(mean.keys()),'values':[y['value'] for y in list(mean.values())]})
        maxval = np.max([np.abs(df_m['values'].min()),np.abs(df_m['values'].max())])
        ymin,ymax = [-maxval,maxval]
        
        fig = px.bar(df_m,x='comps',y='values',
                         labels={'comps':"",'values':'\u0394 ppm'},
                         color='values',
                         color_continuous_scale='tropic_r',
                         hover_data={'values':False,
                                 'comps':False,
                                 'specie':df_m['comps'],
                                 'PPM':(':.1f',df_m['values'])},
                     #width=500,height=500
                    )
        fig.update_layout(
            plot_bgcolor='rgba(0,0,0,0)',
            paper_bgcolor='rgba(0,0,0,0)',
            hovermode="closest",
            hoverlabel=dict(font_size=16),
            coloraxis_showscale=False
        )
        fig.update_xaxes(showgrid=False,tickangle=-60,tickmode='linear')
        dtick=int(int(ymax-ymin)/10)
        fig.update_yaxes(showgrid=True,tickmode='linear',range=[ymin-2,ymax+2],dtick=dtick)

        resultsgraph = dcc.Graph(figure = fig,
                      animate=True,
                      config={'scrollZoom':True},
                      style={'height':'60rem','width':'100%'}
                        )
        
        return([metadata_table,stats_table,paths_table,init_concs_table,final_concs_table,resultsgraph,None])
        
    

    
            
              

if __name__ == '__main__':
    #app.run_server()
    app.run_server(debug = False)  





