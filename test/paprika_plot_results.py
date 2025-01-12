import json
import matplotlib.pyplot as plt
import numpy as np
import math
import os
import sys

sim_type = sys.argv[1]

def loadJson(j_file):
    with open(j_file, 'r') as f:
        d = json.load(f)
    return d

matrices = loadJson("analysis_{0}.json".format(sim_type))

fig, axs = plt.subplots(nrows=1, ncols=2, sharey=True)

attach_fe = matrices['fraction_fe_matrices']['attach']['1.0'][0]
attach_sem = matrices['fraction_sem_matrices']['attach']['1.0'][0]
attach_len = len(attach_fe)

x = np.array(range(0,attach_len))
y = np.array(attach_fe)
yerr = np.array(attach_sem)
ax = axs[0]
ax.errorbar(x,y,yerr,fmt='o')
ax.set_title('Attach')
ax.set_xlabel("Window Number")
ax.set_ylabel("Free Energy")

attach = y[-1]
attach_err = yerr[-1]
print("attach error:", attach_err)

pull_fe = matrices['fraction_fe_matrices']['pull']['1.0'][0]
pull_sem = matrices['fraction_sem_matrices']['pull']['1.0'][0]
pull_len = len(pull_fe)
x = np.array(range(0,pull_len))
y = np.array(pull_fe)
yerr = np.array(pull_sem)
ax = axs[1]
ax.errorbar(x,y,yerr,fmt='o')
ax.set_title('Pull')
ax.set_xlabel("Window Number")
ax.set_ylabel("Free Energy")

pull = y[-1]
pull_err = yerr[-1]
print("pull error:", pull_err)

error = math.sqrt(attach_err**2+pull_err**2)
print("ref_state_work:", matrices['ref_state_work'])
print("attach:", attach, "pull:", pull)
print("Total Free Energy:", attach+pull+matrices['ref_state_work'], "+-", error)

plt.show()

