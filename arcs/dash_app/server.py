from multiprocessing import Condition 
import dash
#import dash_core_components as dcc
#import dash_html_components as html
import dash_bootstrap_components as dbc
import webbrowser
from threading import Timer
from dash_bootstrap_templates import load_figure_template
import dash_loading_spinners as dls
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
from arcs.analysis import AnalyseSampling
from arcs.traversal import Traversal
import pickle
import warnings
import plotly.figure_factory as ff
import networkx as nx 
import plotly.graph_objects as go
from arcs.dash_app.domino import terminate_when_parent_process_dies

####NOTES
#
#need to make settings "off canvas"
#
def start_dash(host:str,port:int,server_is_started:Condition,file_locations='./'):
    terminate_when_parent_process_dies()
    external_stylesheets = [
         dbc.themes.QUARTZ
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
    
    def _markdown_compound(_string):
        md = []
        for i in _string:
            try:
                int(i)
                md.append('<sub>{}</sub>'.format(int(i)))
            except:
                md.append(i)
        return(''.join(md))
    
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
    g = pickle.load(open(file_locations+'SCAN_graph_temp.p','rb'))
    t = Traversal(graph=g,reactions=file_locations+'SCAN_reactions_temp.p')
    
    graph = dbc.Alert('No Data',color='light')# None #html.P('None')
    table = dbc.Alert('No Data',color='light')# None #html.P('None')
    table2 =dbc.Alert('No Data',color='light')# None #html.P('None')
    table3 =dbc.Alert('No Data',color='light')# None #html.P('None')
    table4 =dbc.Alert('No Data',color='light')# None #html.P('None')
    table5 =dbc.Alert('No Data',color='light')# None #html.P('None')
    meta = dbc.Alert('No Data',color='light')#None #html.P('None')
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
    
    
    backgroundcolours='rgba(100,100,120,0.5)'
    
    boxstyle = {'align':'flex',
                #'boxShadow': '1px 1px 1px 1px', 
                #'border-radius': '5px',
                'padding':'2%',
                'backgroundColor':backgroundcolours }
               #}
    
    width={'offset':1,'size':5}
    width={'size':5}
    
    rowclass = 'Row'
    
    ######################layout of DASH template########################
    app = dash.Dash(__name__,external_stylesheets=external_stylesheets)
    #we want 3 columns 
    #column 1 = input concentrations
    #column 2 = graph 
    #column 3 = table of concentration changes

    # in the overlay we want the settings as well as the system information 
    # ultimately a print to file option as a report 

    # in row 2 we will just have the frequent reactions, as well as the paths, (and possibly a click function to see the gibbs free energies of reactions)
    app.layout = html.Div(
        style={'padding':'2rem'},
        children=[
#        dbc.Spinner(
#            id="loading-1", 
#            color='rgba(166, 38, 68,1)',
#            type='border',
#            fullscreen=True,
#            children=html.Div(id="loading-output-1"),
#            fullscreen_style={'background-color':'rgba(0.1,0.1,0.1,0.1)'}
#            ),
        dls.Triangle(
            id="loading-1", 
            color='rgba(166, 38, 68,1)',
            fullscreen=True,
            children=html.Div(id="loading-output-1"),
            fullscreen_style={'background-color':'rgba(0.1,0.1,0.1,0.2)'}
            ),
        dbc.Row(
            className=rowclass,
            style={'margin':'1rem','display':'flex','justify-content': 'space-between','flex-wrap': 'wrap'},
            children=[
                ###### Row 1 Col 1 
                dbc.Col(style=boxstyle,
                         children=[
                             html.H5('Control'),
                             dbc.Row(
                                 align='start',
                                 className=rowclass,
                                 children=sliders
                             ),
                             dbc.Row(
                                 className=rowclass,
                                 align='center',
                                 style={'padding':'2rem',
                                       'align':'center'},
                                 justify=True,
                                 children = [
                                     dbc.Col(
                                         children = [
                                             html.H6('Input Concentrations (PPM)'),
                                             html.Div(
                                                 style={'padding':'0.5rem'},
                                                 children=dcc.Textarea(id='concs-out',
                                                          value="SO2=10,\nNO2=50,\nH2S=30,\nH2O=20",
                                                          style={
                                                              'backgroundColor':backgroundcolours,
                                                              'width':'10rem',
                                                              'height':'21rem',
                                                              'color':'white',
                                                          'border':'2px white solid'},
                                                                       className='form-control')),
                                             dbc.Button(
                                                 style={'border':'2px white solid',
                                                        'marginTop':20,
                                                          'marginLeft':10},
                                                        children='update',
                                                        id='concs-submit',
                                                        n_clicks=0,
                                                        className="me-1",
                                                        color='Primary'),
                                         ]
                                     ),
                                     #########################settings
                                     dbc.Col(
                                         style={'align':'center'},
                                         children = [
                                             html.H6('Algorithm Settings'),
                                             html.Div(
                                                 dbc.Accordion(
                                                     start_collapsed=True,
                                                     children=[
                                                         dbc.AccordionItem(
                                                             ["Number of samples used to get a stochastic median average. Recommended amount >100",
                                                                 dcc.Input(id='samples',
                                                                           value="10",
                                                                           debounce=True,
                                                                           style={'backgroundColor':backgroundcolours,
                                                                                  'color':'white'
                                                                                  },
                                                                           className='form-label mt-4'),
                                                                 ],
                                                                 title='Samples',
                                                                 className='accordion',
                                                             ),
                                                       dbc.AccordionItem(
                                                           ["The maximum path depth is the upper limit for how far down a reaction network a sample goes. Recommended amount ~5",
                                                                     dcc.Input(id='pathdepth',
                                                                               value="5",
                                                                               debounce=True,
                                                                               style={'backgroundColor':backgroundcolours,
                                                                                      'color':'white',
                                                                                      },
                                                                               className='form-label mt-4'
                                                                              )
                                                                     ],
                                                           title='Maximum Path Depth',
                                                           className='accordion'
                                                           ),
                                                       dbc.AccordionItem(
                                                           children=["The discovery cutoff, or probability cutoff, determines the amount at which a certain substance is used in the reactions relative to the largest constituents. This is to weight the reactions in terms of the larger concentrations dominating the reaction probabilities. Recommended amount ~3%",
                                                                     dcc.Input(id='probability_cutoff',
                                                                               value="3",
                                                                               debounce=True,
                                                                               style={'backgroundColor':backgroundcolours,
                                                                                      'color':'white',
                                                                                      },
                                                                               className='form-label mt-4'
                                                                               )
                                                                     ],
                                                           title='Discovery % Cutoff',
                                                           className='accordion'
                                                           ),
                                                       dbc.AccordionItem(
                                                           children=[
                                                                     dcc.Input(id='ceiling',
                                                                               value="2000",
                                                                               debounce=True,
                                                                               style={'backgroundColor':backgroundcolours,
                                                                                      'color':'white'},
                                                                               className='form-label mt-4',
                                                                              ),
                                                           ],
                                                           title='Concentration % Ceiling',
                                                           className='accordion'
                                                           ), 
                                                       dbc.AccordionItem(
                                                           children=[
                                                                     dcc.Input(id='scale_highest',
                                                                               value="0.1",
                                                                               debounce=True,
                                                                               style={
                                                                                      'backgroundColor':backgroundcolours,
                                                                                      'color':'white',
                                                                                      },
                                                                               className='form-label mt-4'
                                                                              )
                                                                     ],
                                                           title='Scale Large Concentrations',
                                                           className='accordion'
                                                           ),
                                                       dbc.AccordionItem(
                                                               children=[
                                                                   dcc.Input(id='max_rank',
                                                                             value="10",
                                                                             debounce=True,
                                                                             style={'backgroundColor':backgroundcolours,
                                                                                    'color':'white',
                                                                                    },
                                                                             className='form-label mt-4'
                                                                             )
                                                                   ]
                                                               ),
                                                       dbc.AccordionItem(
                                                               children=[
                                                                   dcc.Input(id='max_compounds',
                                                                             value="5",
                                                                             debounce=True,
                                                                             style={'backgroundColor':backgroundcolours,
                                                                                    'color':'white',
                                                                                    },
                                                                             className='form-label mt-4'
                                                                             )
                                                                   ]
                                                               ),    
                                                       dbc.AccordionItem(
                                                               children=[
                                                                   dcc.Dropdown(
                                                                       id='method',
                                                                       options={'Bellman-Ford': 'Bellman-Ford',
                                                                                'Dijkstra': 'Dijkstra'},
                                                                       value='Bellman-Ford',
                                                                       className='form-label mt-4')
                                                                       ]
                                                                   ),
                                                       dbc.AccordionItem(
                                                                   children=[
                                                                       dcc.Dropdown(
                                                                           id='include_co2',
                                                                           options={'True':'True',
                                                                                    'False':'False'},
                                                                           value='False',
                                                                           className="form-label mt-4")
                                                                           ]
                                                                        ),
                                                                ],
                                                    ),
                                             ),
                                             html.Div(
                                                 children=dbc.Button(style={'border':'2px white solid',
                                                                            'marginTop':30},
                                                                     children='Submit',
                                                                     id='submit-val',
                                                                     n_clicks=0,
                                                                     className="btn btn-warning",
                                                                     outline=False,
                                                                     color='warning')
                                             ),
                                     ],
                             ),
                             ############################ Row 2
                             dbc.Row(
                                 align='end',
                                 children=[
                                    dbc.Col(children=[
                                        dbc.Card(
                                            dbc.CardBody([
                                                 html.H5('System Information'),
                                                 html.Div(
                                                     id = 'metadata',
                                                     style={'align':'center'},
                                                     children=table3,
                                                     className='table table-dark'
                                                      ),
                                                 ],
                                                         className='card text-white bg-dark mb-3'
                                                         )
                                            )

                                     ],
                                             width={'offset':0,'size':5},
                                            ),
                                     dbc.Col(
                                         children=[
                                             dbc.Card(
                                                 dbc.CardBody([
                                                     html.H5('Concentration Evolution'),
                                                     html.Div(
                                                         id = 'final_concs',
                                                         style={'align':'center'},
                                                         children=table5,
                                                         className='table table-success'
                                                         )],
                                                              className='card text-black bg-success mb-3'
                                                              )
                                                 )
                                         ],
                                         width={'offset':1,'size':6}
                                     )
                                 ],
                             ),
                         ],
                       ),
                #2nd column
                dbc.Col(
                        style=boxstyle,
                        children=[
                            dbc.Card(
                                dbc.CardBody(
                                    [
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
                                    ),
                                ),
                            ],
                        width=width,
                        ),
            ],
        ),
            dbc.Row(
                className=rowclass,
                style={'margin':'1rem','display':'flex','justify-content': 'space-between','flex-wrap': 'wrap'},
                justify=True,
                children=[
                    dbc.Col(style=boxstyle,
                        children=[
                                dbc.Card(
                                    dbc.CardBody([html.H5('Most Frequent Reactions'),
                                                   html.Div(
                                                       id = 'reaction-stats',
                                                       style={'align':'center'},
                                                       children=table,
                                                       className='table table-primary')
                                                   ],
                                                 className='card text-black bg-primary mb-3'
                                                  )
                                    )
                        ],
                            width={'size':6},
                           ),#width={'offset':'1rem'}),
                    dbc.Col(style=boxstyle,
                        children=[
                            dbc.Card(
                                dbc.CardBody([html.H5('Most Frequent Paths'),
                                              html.Div(
                                                  id = 'reaction-paths',
                                                  style={'align':'center'},
                                                  children=table2,
                                                  className='table text-black table-primary'
                                                  ),
                                              ],
                                              className='card text-black bg-primary mb-3'
                                             ),
                                ),
                            ],
                        width={'size':6}
                            ),
                    ],
                    ),
            dbc.Row(
                className=rowclass,
                justify=True,
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
        
    ]),
        ],
    )
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
        settings['probability_threshold'] = float(inputs[2])/100
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
                         'overflowY': 'auto'},
                fixed_rows={'headers': True},
                style_cell_conditional=[
                        {'if': {'column_id': 'index'},
                         'width': '50%',
                        'textAlign':'left'},
                        {'if': {'column_id': ['k','value']},
                         'width': '50%',
                        'textAlign':'left'}],
                                        )
    
            df_d = pd.DataFrame(t.initfinaldiff[ambient_settings['T']][ambient_settings['P']]).round(1).drop('CO2')
            df_d = df_d.T
            df_d  = df_d.T.rename({x:_markdown_compound(x) for x in df_d})
            df_d = df_d.reset_index()
            df_d = df_d.rename({'index':'compound','initial':'initial (ppm)','final':'final (ppm)','change':'change (ppm)'})
    
            
            diff_table = dash_table.DataTable(
                columns=[{'name':i,'id':i,'type':'text','presentation':'markdown'} for i in df_d.columns],
                data=df_d.to_dict('records'),
                style_as_list_view=True,
                cell_selectable=False,
                style_cell={'font_family':'helvetica',
                            'align':'center'},#       'hover-background-color': '#555555'},
                #style_header={'backgroundColor': 'Transparent'},
                #style_data={'backgroundColor': 'Transparent',
                #            'border':'none',
                #           'whiteSpace':'normal'},
                style_table={'height': '400px', 
                         'overflowY': 'auto',
                             },#            'width':'50%'},
                fixed_rows={'headers': True},
                style_cell_conditional=[
                        {'if': {'column_id': 'index'},
                         'width': '10%',
                        'textAlign':'left'},
                        {'if': {'column_id': ['initial','final','change']},
                         'width': '10%',
                        'textAlign':'left'}],
                markdown_options={"html": True, "link_target": "_self"},
                
                                        )
            
            analyse =  AnalyseSampling(t.data,markdown=True)
            analyse.mean_sampling() ; mean = analyse.mean_data[ambient_settings['T']][ambient_settings['P']]
            analyse.reaction_paths() ; paths = pd.DataFrame(analyse.common_paths[ambient_settings['T']][ambient_settings['P']])
            analyse.reaction_statistics() ; stats = pd.DataFrame(analyse.stats[ambient_settings['T']][ambient_settings['P']])
    
            stats_table = dash_table.DataTable(
                columns = [{'name':'Reactions','id':'index','type':'text','presentation':'markdown'},
                               {'name':'k','id':'k','type':'text','presentation':'markdown'},
                               {'name':'Frequency','id':'frequency','type':'text','presentation':'markdown'}],
                data=stats.to_dict('records'),
                style_as_list_view=True,
                cell_selectable=False,
                style_cell={'font_family':'helvetica',
                            'text_align':'center'},
                style_header={'backgroundColor': 'Transparent'},
                style_data={'backgroundColor':'Transparent',
                            'border':'none',
                            'whiteSpace': 'normal'},
                style_table={'align':'centre',
                             'height': '40rem',
                             'width':'100%',
                            'overflowY':'auto'},
                fixed_rows={'headers': True},
                style_cell_conditional=[
                        {'if': {'column_id': 'index'},
                         'width': '60%',
                        'textAlign':'left'},
                        {'if': {'column_id': ['k','frequency']},
                         'width': '20%',
                        'textAlign':'left'}],
                markdown_options={"html": True, "link_target": "_self"},
            )
    #        stats_table = dbc.Table.from_dataframe(stats, striped=True,bordered=False, hover=True, index=False,
    #                                               style={'height':'40rem',
    #                                                      'width':'80%',
    #                                                      'overflow':'auto'},
    #                                               className='overflow-wrapper') 
            paths_table = dash_table.DataTable(
                    columns = [{'name':'Paths','id':'paths','type':'text','presentation':'markdown'},
                               {'name':'k','id':'k','type':'text','presentation':'markdown'},
                               {'name':'Frequency','id':'frequency','type':'text','presentation':'markdown'}],
            #columns=[{'name':i,'id':i} for i in new_stats.columns],
                data=paths.to_dict('records'),
                style_as_list_view=True,
                cell_selectable=False,
                style_cell={'font_family':'helvetica',
                            'text_align':'center',
                            'whiteSpace': 'pre-line',
                            'padding': '3px'},
                style_header={'backgroundColor': 'Transparent'},
                style_data={'backgroundColor': 'Transparent',
                            'border':'none',
                           'whiteSpace': 'break-spaces'},
                style_table={'height': '400px', 
                             'overflowY': 'auto'},
                fixed_rows={'headers': True},
                style_cell_conditional=[
                        {'if': {'column_id': 'Paths'},
                         'width': '30%',
                        'textAlign':'left'},
                        {'if': {'column_id': ['k']},
                         'width': '20%',
                        'textAlign':'left'},
                        {'if': {'column_id': ['Frequency']},
                         'width': '20%',
                        'textAlign':'left'}],
                markdown_options={'html':True,"link_target":"_self"}
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
            
            return([metadata_table,stats_table,paths_table,diff_table,resultsgraph,None])
            
    
    with server_is_started:
        server_is_started.notify()
    
    app.run_server(debug=False,port=port,host=host)  
