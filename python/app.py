import dash
from dash.dependencies import Input, Output
from dash_html_components.Div import Div
import dash_table
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import pandas as pd
import pathlib
import datetime


fname = pathlib.Path('details.csv')
mtime = datetime.datetime.fromtimestamp(fname.stat().st_mtime)
last_modfied_time = mtime.strftime("%d-%b-%Y, %H:%M:%S")

df1 = pd.read_csv('details.csv',sep="\t")
df2 = pd.read_csv('lineage_data.tsv',sep="\t")
df=df1.merge(df2,on='lineage')[['Primer','Direction','lineage','dg_sequence','description','countries']]
df.columns=['Primer','Direction','Lineage','DeltaG','Description','Countries']
#app = dash.Dash(__name__)
last_updated=f"#last Updated {last_modfied_time} based on {len(df['Lineage'].unique())} variants and {len(df['Primer'].unique())} primer sets."



app = dash.Dash(__name__,external_stylesheets=[dbc.themes.BOOTSTRAP])
app.index_string = """<!DOCTYPE html>
<html>
    <head>
        <!-- Global site tag (gtag.js) - Google Analytics -->
        <script async src="https://www.googletagmanager.com/gtag/js?id=G-WLCC056XKX"></script>
        <script>
          window.dataLayer = window.dataLayer || [];
          function gtag(){dataLayer.push(arguments);}
          gtag('js', new Date());

          gtag('config', 'G-WLCC056XKX');
        </script>
        {%metas%}
        <title>{%title%}</title>
        {%favicon%}
        {%css%}
    </head>
    <body>
        {%app_entry%}
        <footer>
            {%config%}
            {%scripts%}
            {%renderer%}
        </footer>
    </body>
</html>"""
app.title = "assayM"

server = app.server

row = html.Div(
    [
        dbc.Row(
            [
                dbc.Col(html.A(href="https://www.biorxiv.org/content/10.1101/2020.12.18.423467v1",children=[html.Img(src=app.get_asset_url('logo1.png'),height=130,width=300)]),width=3),
                dbc.Col(html.Div([html.H1("assayM - track sensitivity of RT-PCR primers on COVID-19 Variants"),
                html.H4(["Developed @ Arnab Pain Lab, enabled by data from",html.A(href="https://www.gisaid.org/",children=[html.Img(src=app.get_asset_url('gisaid.png'),height=30,width=85)]),"and cov-lineages.org"]),
                html.H5(last_updated),html.Div(html.H6("Highlighted DeltaG > 0 implies potential loss of sensitivity of the PCR primer for that specific variant"),style={'color': 'tomato', 'fontSize': 14})]),
                width=9),
            ]
        ),   
    ]
)

ack=html.Div(
[
    dbc.Row(
        [
            dbc.Col(dcc.Markdown('''
            GISAID data provided on this website are subject to GISAID’s [Terms and Conditions](https://www.gisaid.org/DAA/)  


            We would like to thank the GISAID Initiative and are grateful to all of the data contributors, i.e. the Authors, the Originating laboratories responsible for obtaining the specimens, and the Submitting laboratories for generating the genetic sequence and metadata and sharing via the GISAID Initiative, on which this research is based.     
                 
            Elbe, S., and Buckland-Merrett, G. (2017) Data, disease and diplomacy: GISAID’s innovative contribution to global health. Global Challenges, 1:33-46. DOI: [10.1002/gch2.1018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6607375/) PMCID: [31565258](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6607375/)     
                 
            Note: When using results from these analyses in your manuscript, ensure that you also acknowledge the Contributors of data, i.e. “We gratefully acknowledge all the Authors, the Originating laboratories responsible for obtaining the specimens, and the Submitting laboratories for generating the genetic sequence and metadata and sharing via the GISAID Initiative, on which this research is based.”     
            Also, cite the following reference:      
            Shu, Y., McCauley, J. (2017) GISAID: From vision to reality. EuroSurveillance, 22(13) DOI: [10.2807/1560-7917.ES.2017.22.13.30494](http://dx.doi.org/10.2807/1560-7917.ES.2017.22.13.30494) PMCID: [PMC5388101](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5388101/)
            '''))
        ]
    ),
]
)


app.layout=html.Div([
    row,
    dash_table.DataTable(
        id='datatable-interactivity',
        columns=[
            {"name": i, "id": i, "deletable": False, "selectable": True} for i in df.columns
        ],
        data=df.to_dict('records'),
        style_cell={
        'whiteSpace': 'normal',
        'height': 'auto',
        'font-family':'sans-serif'
        },
        style_header_conditional=[
        {'if': {'column_id': 'Primer'},
         'width': '10%'},
        {'if': {'column_id': 'Lineage'},
         'width': '10%'},
            {
        'if': {'column_id': 'Primer'},
        'backgroundColor': 'darkturquoise',
        'color': 'black'
        },
        {
        'if': {'column_id': 'Direction'},
        'backgroundColor': 'darkturquoise',
        'color': 'black'
        },
        {
        'if': {'column_id': 'Lineage'},
        'backgroundColor': 'olive',
        'color': 'black'
        },
        {
        'if': {'column_id': 'DeltaG'},
        'backgroundColor': 'tomato',
        'color': 'white'
        },
        {
        'if': {'column_id': 'Description'},
        'backgroundColor': 'deepskyblue',
        'color': 'black'
        },
                {
        'if': {'column_id': 'Countries'},
        'backgroundColor': 'darkslateblue',
        'color': 'white'
        }
        ],
        style_data_conditional=[{
            'if': {
                'filter_query': '{DeltaG} > 0',
                'column_id': ['Primer','Direction','Lineage','DeltaG','Description','Countries'],
            },
            'backgroundColor': 'tomato',
            'color': 'white'
        }],
        editable=False,
        filter_action="native",
        sort_action="native",
        sort_mode="multi",
        #column_selectable="single",
        #row_selectable="multi",
        row_deletable=False,
        selected_columns=[],
        selected_rows=[],
        page_action="native",
        page_current= 0,
        page_size= 50,
    ),
    html.Div(id='datatable-interactivity-container'),
    ack
])




if __name__ == '__main__':
    app.run_server(debug=False)
