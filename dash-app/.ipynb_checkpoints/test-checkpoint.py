import dash
#import dash_core_components as dcc
#import dash_html_components as html
import dash_bootstrap_components as dbc
from dash import html
from dash import dcc
import plotly.express as px
import pandas as pd


external_stylesheets = [
     #dbc.themes.YETI
     #dbc.themes.QUARTZ #- funky
     dbc.themes.MINTY
]


app = dash.Dash(__name__,external_stylesheets=external_stylesheets)

df = pd.DataFrame({
  "Indicator": ["Total Cases", "Recovered", "Total Cases", "Recovered"],
  "Cases": [19111326, 11219123, 10147468, 9717834],
  "Country": ["USA", "USA", "India", "India"]
})

fig = px.bar(df, x = "Indicator", y = "Cases", color = "Country", barmode = "group",
             width=800,height=800)

fig.update_layout(
    margin=dict(l=20, r=20, t=20, b=20),
    paper_bgcolor="rgba(0,0,0,0)",
    plot_bgcolor='rgba(0,0,0,0)'
)

app.layout = html.Div(
        style = {'plot_bgcolor':'rgba(0,0,0,0)','font_color':'#FFFFFF'}, 
        children = [
            #first row = header
            html.Div(style={'marginBottom':50,'marginTop':50,'display':'inline-block'},#align vertically
                     children=[
                         html.Img(src=app.get_asset_url('ARCS_Logo-01.png'),
                                 style={'height':'10%', 'width':'10%'},
                                  className = 'center'),
                        # html.H1(children = 'ARCS',
                        #         style = {'textAlign': 'center','align':'center','className':'h-50'}
                        #         ),
                         html.Div(children = ['Automated Reactions for CO',html.Sub(2),' Storage'], 
                                  style = {'textAlign': 'center'}
                                  ),
                         ],className='row'),
            #second row = plot field and sliders
            html.Div(style={'display':'inline-block'},
                            children=[
                                #slider first column
                                html.Div(children=dcc.Slider(-5, 10, 1, value=-3)
                             ),
                                #graph second column
                                dcc.Graph(id = 'example-graph-2',
                                          figure = fig,
                                          config={'scrollZoom':True}
                                          )
                                ],className='row')
                     ])

if __name__ == '__main__':
  app.run_server(debug = True)
