f = open('build/bin/log.txt')
content = f.read()
f.close()

import re
idx = [m.start() for m in re.finditer('HybridInterpolatedCurve', content)]
subs = []

for n in range(0,len(idx)-1):
    subs.append(content[idx[n]:idx[n+1]])

# f = open('values.txt', 'w')
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
            print('skip')
            continue
    except IndexError:
        print('passing')
        pass
    # f.write(t + ' | ' + val + '\n')
    whole += val + '\n'
# f.close()

import numpy as np
from io import StringIO
import scipy.io as sio
values = np.loadtxt(StringIO(whole), delimiter=',')
print(values.shape)

# CHECK THE ORDER
labels = ['q_dot_ref','q_ref','qs_dot','qs','qm_dot','qm','qm_m2s','qm_dot_m2s','qs_s2m','qs_dot_s2m','H_in_m','H_in_s','qm_m2s_p','qm_dot_m2s_p','qs_s2m_p','qs_dot_s2m_p','H_in_m_p',
            'H_in_s_p','qm_d','qs_d','qm_dot_d','qs_dot_d','fe_d','cnt','t','qs_d_prev','qm_d_prev','tau_plm','Hm','H-m','tau_tlc','deltaH_m','tau_pls','Hs','H-s','deltaH_s','x_rand',
            'y_rand','z_rand','h_m_star','tau_tlm','tau_tls','fe']

assert(len(labels) == values.shape[1])

matfile = {}
for val, label in zip(values.T, labels):
    matfile[label]=val
sio.savemat('build/bin/values.mat', matfile)