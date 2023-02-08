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

external_stylesheets = [
     #dbc.themes.YETI
     dbc.themes.QUARTZ, #- funky
     #dbc.themes.MINTY
     #dbc.themes.SUPERHERO
]
#load_figure_template('MINTY')
load_figure_template('QUARTZ')
#load_figure_template('SUPERHERO')


backgroundcolours='rgba(0,0,0,0.1)'

app = dash.Dash(__name__,external_stylesheets=external_stylesheets)

df = pd.DataFrame({
  "Indicator": ["Total Cases", "Recovered", "Total Cases", "Recovered"],
  "Cases": [19111326, 11219123, 10147468, 9717834],
  "Country": ["USA", "USA", "India", "India"]
})

fig = px.bar(df, x = "Indicator", y = "Cases", color = "Country", barmode = "group",
             width=500,height=500)

fig.update_layout(
    #margin=dict(l=50, r=50, t=50, b=50),
    plot_bgcolor='rgba(0,0,0,0)',
    paper_bgcolor='rgba(0,0,0,0)',
    #xaxis={'showline':True,'linewidth':2,'mirror':True},
    #yaxis={'showline':True,'linewidth':2,'mirror':True}
)

app.layout = html.Div(
    style={'padding':'4rem'},
    children=[
    
    #first row = header
    dbc.Row(
        style={'padding':'1rem',
               'marginBottom':'1%',
               'marginTop':'1%',
               'display':'flex',
               'width':'100%',
               'align':'center'},
        children=[          
            dbc.Col(
                html.Div(style={'display':'flex','align':'end'},
                             children=[
                                 html.Img(src=app.get_asset_url('ARCS_Logo-01.png'),
                                          style={'width':'10%','height':'10%','padding':'1rem','align':'end'}),
                                          #className = 'center'),
                             ]
                        )),
        ],className='row1-headers'),
    #second row = plot field and sliders
    dbc.Row(
        style={'margin':'1rem','display':'flex','justify-content': 'space-between','flex-wrap': 'wrap'},
        children=[
            #1st column
            dbc.Col(style={'align':'center',
                           'boxShadow': '2px 2px 2px 2px', 
                           'border-radius': '10px',
                           'padding':'1rem', 
                           'backgroundColor':backgroundcolours},
                    children=[dcc.Graph(id = 'graph',
                              figure = fig,
                                        animate=True,
                              config={'scrollZoom':False}
                             )],
                   ),
            #2nd column
            dbc.Col(style={'align':'center',
                            'boxShadow': '2px 2px 2px 2px', 
                            'border-radius': '10px',
                            'padding':'2rem',
                            'backgroundColor':backgroundcolours},
                     children=[
                         html.Div(
                             #style={'columnCount':2},
                             children=[
                                 html.Label(style={'width':'100%','textAlign':'left','paddingLeft': '1rem'},
                                        children=['CO',html.Sub(2)]),
                                 html.Div(style={'width':'90%'},
                                          children=dcc.Slider(-5, 10, 1, value=0,
                                                              id='slider1'),
                                          className='center')
                             ]),
                         html.Div(
                             #style={'display':'block'},
                             children=[
                                 html.Label(style={'width':'100%','textAlign':'left','paddingLeft': '1rem'},
                                        children=['SO',html.Sub(2)]),
                                 html.Div(style={'width':'90%'},
                                          children=dcc.Slider(-5, 10, 1, value=0,
                                                              id='slider2'),
                                          className='center')
                             ]),
                         html.Div(
                             #style={'display':'block'},
                             children=[
                                 html.Label(style={'width':'100%','textAlign':'left','paddingLeft': '1rem'},
                                        children=['NO',html.Sub(2)]),
                                 html.Div(style={'width':'90%'},
                                          children=dcc.Slider(-5, 10, 1, value=0,
                                                              id='slider3'),
                                          className='center')
                             ]),
                         html.Div(
                             #style={'display':'block'},
                             children=[
                                 html.Label(style={'width':'100%','textAlign':'left','paddingLeft': '1rem'},
                                        children='Concentrations'),
                                 html.Div(style={'width':'90%',
                                                 'display':'flex',
                                                 'justifyContent':'center',
                                                 'color':'white'},
                                          children=[dcc.Textarea(id='input1',
                                                             value="CO2:100,\nSO2:100,\nNO2:50",
                                                             style={'marginRight':'10px',
                                                                    'backgroundColor':backgroundcolours,
                                                                    'width':'100%',
                                                                    'height':200,
                                                                   'color':'white'}),
                                                    dbc.Button(style={'height':'10%'},
                                                                children='Submit',
                                                                id='submit-val',
                                                               n_clicks=0,
                                                               className="me-1",
                                                               color='Secondary')],
                                          className='center'),
                             ]),
                     ],
                   width={'offset':1})
        ],className='row'),
        dbc.Row(
            children=[html.Img(src=app.get_asset_url('ntnu-inv.png'),
                                          style={'width':'15%','height':'15%','padding':'2rem','align':'center'}),
                     ]
        )
    
])

#@app.callback(
#    Output("output", "children"),
#    Input("submit-val", "n_clicks"),
#    State("input1","value")
#)
#def update_output(input1):
#    if input1 == None:
#        concs = {'CO2':100,'SO2':100,'NO2':50}
#    concs = {}
#    itext = str(input1).split(',')
#    for i in itext:
#        spec,num = i.split('=')
#        concs[spec] = float(num)
#    return(pd.Series(concs))

if __name__ == '__main__':
  app.run_server(debug = True)

