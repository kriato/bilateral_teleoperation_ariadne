f = open('build/bin/log.txt')
content = f.read()
f.close()

import re
idx = [m.start() for m in re.finditer('HybridInterpolatedCurve', content)]
subs = []

for n in range(0,len(idx)-1):
    subs.append(content[idx[n]:idx[n+1]])

# f = open('values.txt', 'w')
whole = ''
for episode in subs:
    points_idx = episode.index('points')
    regex = re.compile(r'\[.*?\]')
    res = re.findall(regex, episode[points_idx:])
    for tmp in res:
        # f.write(tmp[1:-1] + '\n')
        whole += tmp[1:-1] + '\n'

# f.close()

import numpy as np
from io import StringIO
import scipy.io as sio
values = np.loadtxt(StringIO(whole), delimiter=',')
print(values.shape)

# CHECK THE ORDER
labels = ['q_dot_ref','q_ref','qs_dot','qs','qm_dot','qm','qm_d','qs_d','qm_dot_d','qs_dot_d','fe_d','cnt','t','qs_d_prev','qm_d_prev','tau_plm','Hm','Hm_out','tau_tlc','deltaH_m','tau_pls','Hs','Hs_out','deltaH_s','fe','qm_m2s','qm_dot_m2s','tau_tls','qs_s2m','qs_dot_s2m','h_m_star','tau_tlm','H+m','H+s']

assert(len(labels) == values.shape[1])

matfile = {}
for val, label in zip(values.T, labels):
    matfile[label]=val
sio.savemat('build/bin/values.mat', matfile)