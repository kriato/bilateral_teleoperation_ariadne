f = open('build/bin/log.txt')
content = f.read()
f.close()

import re
idx = [m.start() for m in re.finditer('HybridInterpolatedCurve', content)]
subs = []

for n in range(0,len(idx)-1):
    subs.append(content[idx[n]:idx[n+1]])

# f = open('values.csv', 'w')
all_ts = []
all_val = []
whole = ''
for episode in subs:
    points_idx = episode.index('points')
    regex = re.compile(r'\[.*?\]')
    res = re.findall(regex, episode[points_idx:])
    regex = re.compile(r' .*?\:')
    ts = re.findall(regex, episode[points_idx:])
    for index, (t,val) in enumerate(zip(ts,res)):
        all_ts.append(t[1:-1])
        all_val.append(val[1:-1])

assert(len(all_ts), len(all_val))
for index, (t,val) in enumerate(zip(all_ts,all_val)):
    try:
        if t == all_ts[index+1]:
            # print('skip')
            continue
    except IndexError:
        # print('passing')
        pass
    # f.write(t + ',' + val + '\n')
    whole += val + '\n'
# f.close()

import numpy as np
from io import StringIO
import scipy.io as sio
values = np.loadtxt(StringIO(whole), delimiter=',')
print(values.shape)

# CHECK THE ORDER
labels = ['q_dot_ref','q_ref','qs_dot','qs','qm_dot','qm','qm_d','qs_d','qm_dot_d','qs_dot_d','fe_d','cnt','t','qs_d_prev','qm_d_prev','tau_plm','Hm','H-m','tau_tlc','deltaH_m','tau_pls','Hs',
            'H-s','deltaH_s','x_rand','y_rand','z_rand','rand','fe','qm_m2s','qm_dot_m2s','tau_tls','qs_s2m','qs_dot_s2m','h_m_star','tau_tlm','H+m','H+s']

assert(len(labels) == values.shape[1])

matfile = {}
for val, label in zip(values.T, labels):
    matfile[label]=val
sio.savemat('build/bin/values.mat', matfile)

print(max(matfile['Hm']))
#Import plotly

import plotly.graph_objs as go
import plotly.express as px
import pandas as pd
df = pd.DataFrame(matfile)

#Initialize the plot
# fig_totalsales=px.
# #Add the lines
# fig_totalsales.add_trace(go.(x=sales[sales.columns[0]], y=sales[sales.columns[1]], visible=True, name='Sales'))
# fig_totalsales.add_trace(go.Scatter(x=fore['Date'],y=fore[fore.columns[0]],visible=True, name='Forecast'))
# #Add title/set graph size
# fig_totalsales.update_layout(title = 'Daily Sales', width=850, height=400)
# fig_totalsales.show()
fig = px.line(df, x='t',  y=['qm','qs','q_ref'])
fig.update_layout(title = {
    'text': '<b>Master/Slave position vs reference</b>'
    # 'xanchor': 'center',
    # 'yanchor': 'top'
    }, 
    legend_title='<b>Legend</b>',
    yaxis_title='rad',
    xaxis_title='Time')

fig_2 = px.line(df, x='t',  y=['qm_dot','qs_dot','q_dot_ref'])
fig_2.update_layout(title = {
    'text': '<b>Master/Slave velocity vs reference</b>'
    # 'xanchor': 'center',
    # 'yanchor': 'top'
    }, 
    legend_title='<b>Legend</b>',
    yaxis_title='rad/s',
    xaxis_title='Time')

fig_3 = px.line(df, x='t',  y=['Hm','Hs'])
fig_3.update_layout(title = {
    'text': '<b>Master and Slave energy tank levels</b>'
    # 'xanchor': 'center',
    # 'yanchor': 'top'
    }, 
    legend_title='<b>Legend</b>',
    yaxis_title='Joule',
    xaxis_title='Time')

fig_4 = px.line_3d(df, x='x_rand', y='y_rand' ,z='z_rand', labels={'x_rand': 'x', 'y_rand': 'y', 'z_rand': 'z'})
fig_4.update_layout(title = {
    'text': '<b>Lorentz system position</b>'
    # 'xanchor': 'center',
    # 'yanchor': 'top'
    }, 
    legend_title='<b>Legend</b>',
    yaxis_title='Meters',
    xaxis_title='Time')

import dash
import dash_core_components as dcc
import dash_html_components as html

app = dash.Dash()
app.layout = html.Div(children=[
    html.H1(children='Bilateral teleoperation simulated with Hybrid Automata'),
    html.Div([
    dcc.Graph(figure=fig, style={'margin-left': 'auto','margin-right': 'auto', 'display': 'block', 'height': 500, 'width': 1200}),
    dcc.Graph(figure=fig_2, style={'margin-left': 'auto','margin-right': 'auto', 'display': 'block', 'height': 500, 'width': 1200}),
    dcc.Graph(figure=fig_3, style={'margin-left': 'auto','margin-right': 'auto', 'display': 'block', 'height': 500, 'width': 1200}),
    dcc.Graph(figure=fig_4, style={'margin-left': 'auto','margin-right': 'auto', 'display': 'block', 'height': 600, 'width': 700})
    ])
], style={'textAlign': 'center'})

app.run_server(host='localhost', debug=True)