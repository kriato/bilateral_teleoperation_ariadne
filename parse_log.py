f = open('build/bin/log.txt')
content = f.read()
f.close()

import re
idx = [m.start() for m in re.finditer('HybridInterpolatedCurve', content)]
subs = []

var_list = []
for n in range(0,len(idx)-1):
    subs.append(content[idx[n]:idx[n+1]])

import numpy as np
# f = open('values.csv', 'w')
all_ts = []
all_val = []
all_vars = []
all_sort_indices = []
whole = ''
for episode in subs:
    regex = re.compile(r'\[.*?\]')
    res = re.findall(regex, episode)
    variables = res[0][1:-1].split(',')
    sort_idx = np.argsort(variables)
    points_idx = episode.index('points')
    regex = re.compile(r'\[.*?\]')
    res = re.findall(regex, episode[points_idx:])
    regex = re.compile(r' .*?\:')
    ts = re.findall(regex, episode[points_idx:])
    for index, (t,val) in enumerate(zip(ts,res)):
        all_ts.append(t[1:-1])
        all_val.append(val[1:-1])
        all_vars.append([variables[i] for i in sort_idx])
        all_sort_indices.append(np.asarray(sort_idx))

final_ordered_vars = []
final_all_sort_indices = []
assert(len(all_ts), len(all_val), len(all_vars), len(all_sort_indices))
for index, (t,val,var,sort_index) in enumerate(zip(all_ts,all_val,all_vars,all_sort_indices)):
    try:
        if t == all_ts[index+1]:
            # print('skip')
            continue
    except IndexError:
        # print('passing')
        pass
    # f.write(t + ',' + val + '\n')
    whole += val + '\n'
    final_ordered_vars.append(var)
    final_all_sort_indices.append(sort_index)
# f.close()

from io import StringIO
import scipy.io as sio
values = np.loadtxt(StringIO(whole), delimiter=',')
print(values.shape)

# CHECK THE ORDER
labels = ['H+m', 'H+s', 'H-m', 'H-s', 'Hm', 'Hs', 'cnt', 'deltaH_m', 'deltaH_s', 'fe', 'fe_d', 'h_m_star', 'q_dot_ref', 'q_ref', 'qm', 'qm_d', 'qm_d_prev', 'qm_dot', 'qm_dot_d', 'qm_dot_m2s', 'qm_m2s', 'qs', 'qs_d', 'qs_d_prev', 'qs_dot', 'qs_dot_d', 'qs_dot_s2m', 'qs_s2m', 't', 'tau_plm', 'tau_pls', 'tau_tlc', 'tau_tlm', 'tau_tls', 'x_rand', 'y_rand', 'z_rand']

for var in final_ordered_vars:
    assert(var, labels)

assert(len(labels) == values.shape[1])

ordered_values = np.zeros(values.shape)
assert(values.shape[0], len(final_all_sort_indices))
for index, (v,idx) in enumerate(zip(values,final_all_sort_indices)):
    ordered_values[index,:] = v[idx]

mat = {}
for val, label in zip(ordered_values.T, labels):
    mat[label]=val

sio.savemat('./build/bin/values.mat', mat)

import plotly.graph_objs as go
import plotly.express as px
import pandas as pd
df = pd.DataFrame(mat)

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

app.run_server(host='192.168.42.13', port='25565', debug=True)