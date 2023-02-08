import dash
#import dash_core_components as dcc
#import dash_html_components as html
import dash_bootstrap_components as dbc
from dash_bootstrap_templates import load_figure_template
from dash import html
from dash import dash_table
from dash import dcc
from dash.dependencies import Input, Output, State
import plotly.express as px
import pandas as pd
import time
import base64
import io
import numpy as np 
from monty.serialization import loadfn
from arcs.setup_functions import GenerateInitialConcentrations

import pickle

external_stylesheets = [
     #dbc.themes.YETI
     dbc.themes.QUARTZ, #- funky
    #dbc.themes.CYBORG
     #dbc.themes.MORPH
     #dbc.themes.QUARTZ
     #dbc.themes.MINTY
     #dbc.themes.SUPERHERO
]
#load_figure_template('YETI')
load_figure_template('QUARTZ')
#load_figure_template('CYBORG')
#load_figure_template('MINTY')
#load_figure_template('MORPH')
#load_figure_template('SUPERHERO')

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
                                className='./assets/slider.css'
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

g = pickle.load(open('./assets/graph.p','rb'))

backgroundcolours='rgba(0,0,0,0.1)'
graph = None
table = None
table2 = None
sliders = make_sliders(g,labels={0:'T',1:'P'})
concs = GenerateInitialConcentrations(g).from_file('./assets/initial_concentrations.json')







######################starting#######################
app = dash.Dash(__name__,external_stylesheets=[external_stylesheets[0],'./assets/table.css'])

app.layout = html.Div(
    style={'padding':'2rem'},
    children=[
    dbc.Row(
        className='Row',
        style={'margin':'1rem','display':'flex','justify-content': 'space-between','flex-wrap': 'wrap'},
        children=[
            #1st column
            dbc.Col(
                    style={'align':'center',
                           'boxShadow': '2px 2px 2px 2px pink', 
                           'border-radius': '5px',
                           'padding':'1rem', 
                           'backgroundColor':backgroundcolours},
                    children=[
                        html.H5('Data'),
                        html.Div(id='output-graph',
                                 children=graph
                                ),
                        html.Div(id='placeholder', 
                                 children=None,
                                 style={'display':'none'})
                             ]
                
                   ),
            #2nd column
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
                                                      value="SO2=100,\nNO2=50,\nH2S=30",
                                                      style={
                                                          'backgroundColor':backgroundcolours,
                                                          'width':'10rem',
                                                          'height':'7rem',
                                                          'color':'white',
                                                      'border':'2px white solid'}))
                                     ]
                                 ),
                                 dbc.Col(
                                     style={'align':'center'},
                                     children = [
                                         html.H6('Algorithm Settings'),
                                         html.Div(
                                             style={'padding':'0.5rem'},
                                             children=dcc.Input(id='samples',
                                                       value="samples",
                                                       style={'backgroundColor':backgroundcolours,
                                                             'color':'white',
                                                             'border':'2px white solid'}
                                                  )),
                                         html.Div(
                                             style={'padding':'0.5rem'},
                                             children=dcc.Input(id='pathdepth',
                                                           value="path depth",
                                                           style={'backgroundColor':backgroundcolours,
                                                                 'color':'white',
                                                                  'border':'2px white solid'}

                                                          )),
                                         html.Div(
                                             style={'padding':'0.5rem'},
                                             children=dcc.Input(id='percentage cutoff',
                                                           value="percentage cutoff",
                                                           style={'backgroundColor':backgroundcolours,
                                                                 'color':'white',
                                                                  'border':'2px white solid'}
  
                                                          )), 
                                         ]
                                 ),
                                 dbc.Col(
                                     dbc.Button(style={'border':'2px white solid',
                                                      'marginTop':50},
                                                children='Submit',
                                                id='submit-val',
                                                n_clicks=0,
                                                className="me-1",
                                                color='Primary')
                                 )
                             ]
                         )
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
            children=[
                dbc.Col(
                    children=[
                        html.Div(style={'display':'inline-block','align':'end'},
                                 children=[
                                     html.Img(src=app.get_asset_url('logos.png'),
                                              style={'width':'50%','height':'50%','padding':'0.5rem','align':'end'}),
                                 ]
                                ),
                    ]
                )

            ]
        )
    
])
 
@app.callback(Output('placeholder','children'),
              Input('concs-out','value'))    

def update_output(input1):
    concs_2 = {}
    itext = str(input1).split(',')
    for i in itext:
        spec,num = i.split('=')
        if '\n' in spec:
            spec = spec.split('\n')[1]
        concs[spec] = float(num)*1e-6
    print(concs)
    #return(concs_2)
    
    
if __name__ == '__main__':
    app.run_server(debug = True)  
    
    
####### slider magic - should be slider amount agnostic

#@app.callback(Output('output-graph','children'),
#              [Input('slider-{}'.format(i),'value') for i in range(len(sliders))])
#def update_graph(*v):
#
#    keys = keys_by_depth(data)    
#    di = {i:sorted(list(keys[i]))[0] for i in keys}
#    
#    nv = []
#    for n in [float(x) for x in v]:
#        if n == 0.00:
#            nv.append('{:.1F}'.format(0.00))
#        else:
#            nv.append('{:.1F}'.format(n))
#        
#    if len(nv) == 1:
#        df1 = data[str(nv[0])]
#    elif len(nv) == 2:
#        df1 = data[str(nv[0])][str(nv[1])]
#    elif len(nv) == 3:
#        df1 = data[str(nv[0])][str(nv[1])][str(nv[2])]
#    elif len(nv) == 4:
#        df1 = data[str(nv[0])][str(nv[1])][str(nv[2])][str(nv[3])]
#    elif len(nv) == 5:
#        df1 = data[str(nv[0])][str(nv[1])][str(nv[2])][str(nv[3])][str(nv[4])]
#    elif len(nv) == 6:
#        df1 = data[str(nv[0])][str(nv[1])][str(nv[2])][str(nv[3])][str(nv[4])][str(nv[5])]
#    elif len(nv) == 7:
#        df1 = data[str(nv[0])][str(nv[1])][str(nv[2])][str(nv[3])][str(nv[4])][str(nv[5])][str(nv[6])]
#             
#    #df1 = data[str(nv[0])][str(nv[1])][str(nv[2])][str(nv[3])][str(nv[4])][str(nv[5])] # needs to be more agnostic
#    df2 = pd.DataFrame({'comps':list(df1.keys()),'values':[y/1e-6 for y in list(df1.values)]})
#    
#    maxval = np.max([np.abs(df2['values'].min()),np.abs(df2['values'].max())])
#    ymin,ymax = [-maxval,maxval]
#
#    # for the slider 
#    i1range = list(data)
#           
#    #figure
#    fig = px.bar(df2,x='comps',y='values',
#                    labels={'comps':"",'values':'\u0394 ppm'},
#                     color='values',
#                     color_continuous_scale='tropic_r',
#                     hover_data={'values':False,
#                                 'comps':False,
#                                 'specie':df2['comps'],
#                                 'PPM':(':.1f',df2['values'])},
#                     #width=700,height=700
#                    )
#    fig.update_layout(
#            plot_bgcolor='rgba(0,0,0,0)',
#            paper_bgcolor='rgba(0,0,0,0)',
#            hovermode="closest",
#            hoverlabel=dict(font_size=16),
#            coloraxis_showscale=False
#        )
#    fig.update_xaxes(showgrid=False,tickangle=-60,tickmode='linear')
#    fig.update_yaxes(showgrid=True,tickmode='linear',range=[ymin-2,ymax+2])
#
#    return(dcc.Graph(figure = fig,
#                      animate=True,
#                      config={'scrollZoom':True},
#                      #style={'height':'100%','width':'100%'}
#                        )
#        ) 

#@app.callback(Output('reaction-stats','children'),
#              [Input('slider-{}'.format(i),'value') for i in range(len(sliders))])
#
#def update_stats(*v):
#    
#    nv = []
#    for n in [float(x) for x in v]:
#        if n == 0.00:
#            nv.append('{:.1F}'.format(0.00))
#        else:
#            nv.append('{:.1F}'.format(n))
#    
#    if len(nv) == 1:
#        df1 = stats[str(nv[0])]
#    elif len(nv) == 2:
#        df1 = stats[str(nv[0])][str(nv[1])]
#    elif len(nv) == 3:
#        df1 = stats[str(nv[0])][str(nv[1])][str(nv[2])]
#    elif len(nv) == 4:
#        df1 = stats[str(nv[0])][str(nv[1])][str(nv[2])][str(nv[3])]
#    elif len(nv) == 5:
#        df1 = stats[str(nv[0])][str(nv[1])][str(nv[2])][str(nv[3])][str(nv[4])]
#    elif len(nv) == 6:
#        df1 = stats[str(nv[0])][str(nv[1])][str(nv[2])][str(nv[3])][str(nv[4])][str(nv[5])]
#    elif len(nv) == 7:
#        df1 = stats[str(nv[0])][str(nv[1])][str(nv[2])][str(nv[3])][str(nv[4])][str(nv[5])][str(nv[6])]
#
#
#    new_stats = pd.DataFrame(df1)
#    if new_stats.empty:
#        new_stats = pd.DataFrame({'index':['-' for x in range(10)],
#                     'k':['-' for x in range(10)],
#                     'frequency':['-' for x in range(10)]})
#
#    return(dash_table.DataTable(
#        columns=[{'name':i,'id':i} for i in new_stats.columns],
#        data=new_stats.to_dict('records'),
#        style_as_list_view=True,
#        cell_selectable=False,
#        style_cell={'font_family':'helvetica',
#                    'text_align':'center',
#                   'hover-background-color': '#555555'},
#        style_header={'backgroundColor': 'Transparent'},
#        style_data={'backgroundColor': 'Transparent',
#                    'border':'none'},
#        style_table={'height': '400px', 
#                     'overflowY': 'auto'}
#                                    )
#          )


#@app.callback(Output('reaction-paths','children'),
#              [Input('slider-{}'.format(i),'value') for i in range(len(sliders))])
#
#def update_paths(*v):
#    
#    nv = []
#    for n in [float(x) for x in v]:
#        if n == 0.00:
#            nv.append('{:.1F}'.format(0.00))
#        else:
#            nv.append('{:.1F}'.format(n))
#
#    if len(nv) == 1:
#        df1 = paths[str(nv[0])]
#    elif len(nv) == 2:
#        df1 = paths[str(nv[0])][str(nv[1])]
#    elif len(nv) == 3:
#        df1 = paths[str(nv[0])][str(nv[1])][str(nv[2])]
#    elif len(nv) == 4:
#        df1 = paths[str(nv[0])][str(nv[1])][str(nv[2])][str(nv[3])]
#    elif len(nv) == 5:
#        df1 = paths[str(nv[0])][str(nv[1])][str(nv[2])][str(nv[3])][str(nv[4])]
#    elif len(nv) == 6:
#        df1 = paths[str(nv[0])][str(nv[1])][str(nv[2])][str(nv[3])][str(nv[4])][str(nv[5])]
#    elif len(nv) == 7:
#        df1 = paths[str(nv[0])][str(nv[1])][str(nv[2])][str(nv[3])][str(nv[4])][str(nv[5])][str(nv[6])]
#
#
#    new_paths = pd.DataFrame(df1)
#    if new_paths.empty:
#        new_paths = pd.DataFrame({'paths':['-' for x in range(10)],
#                                  'k':['-' for x in range(10)],
#                                  'frequency':['-' for x in range(10)]})
#        
#    return(dash_table.DataTable(
#        columns = [{'name':'paths','id':'paths'},{'name':'k','id':'k'},{'name':'frequency','id':'frequency'}],
#        #columns=[{'name':i,'id':i} for i in new_stats.columns],
#        data=new_paths.to_dict('records'),
#        style_as_list_view=True,
#        cell_selectable=False,
#        style_cell={'font_family':'helvetica',
#                    'text_align':'center',
#                    'whiteSpace': 'pre-line',
#                   'padding': '5px'},
#        style_header={'backgroundColor': 'Transparent'},
#        style_data={'backgroundColor': 'Transparent',
#                    'border':'none'},
#        style_table={'height': '400px', 
#                     'overflowY': 'auto'}
#                                    )
#          )




