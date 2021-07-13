#---------Import---------
import sys
import time
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
sys.path.append('../lib/')
from models import *
from networks import *
from sim_loops import *
from utilities import *
from Epi_models import *

print('Hola mundo \n')

numNodes = 20
baseGraph    = networkx.barabasi_albert_graph(n=numNodes, m=9)
G_normal     = custom_exponential_graph(baseGraph, scale=100)
# Social distancing interactions:
G_distancing = custom_exponential_graph(baseGraph, scale=10)
# Quarantine interactions:
G_quarantine = custom_exponential_graph(baseGraph, scale=5)

model = ExtSEIRSNetworkModel(G=G_normal, beta=0.4, sigma=10000, lamda = 1/4, gamma = 1/6, gamma_asym=1/6, a = 1,  p=1,
                          G_Q=G_normal, initI_asym=1, isolation_time = 0, sigma_Q = 1000, lamda_Q = 1/4, gamma_Q_asym = 1/6,
                          transition_mode = 'time_in_state', prevalence_ext = 1e-2, eta = 1, mu_H = 1, xi = 1e-10, nu = 1e-10)

checkpoints = {'t': [20, 100], 'G': [G_distancing, G_normal], 'p': [0.1, 0.5], 'theta_E': [0.02, 0.02], 'theta_I': [0.02, 0.02], 'phi_E':   [0.2, 0.2], 'phi_I':   [0.2, 0.2]}

#model.run(T=7*4, checkpoints=None)
run_tti_sim(model = model, T=7*4, testing_cadence = 'semiweekly2', testing_compliance_random = np.array([True]*model.numNodes), isolation_compliance_positive_individual = [True]*model.numNodes,
	isolation_compliance_positive_groupmate = np.array([True]*model.numNodes), isolation_groups = np.array([range(numNodes)]*model.numNodes), 
	isolation_lag_positive = 0)

numQ = model.numQ_S + model.numQ_E + model.numQ_pre + model.numQ_asym + model.numQ_R

fig, ax = plt.subplots(figsize = (12,8))
ax.plot(model.tseries, model.numS, label = 'S')
ax.plot(model.tseries, model.numE, label = 'E')
ax.plot(model.tseries, model.numI_pre, label = 'pre')
ax.plot(model.tseries, model.numI_sym, label = 'I')
ax.plot(model.tseries, model.numI_asym, label = 'I_a')
ax.plot(model.tseries, numQ , label = 'Q_a')
ax.plot(model.tseries, model.numR, label = 'R')
my_plot_layout(ax = ax)
ax.legend()
plt.show()

fig2, ax2 = plt.subplots(figsize = (12,8))
#ax2.stackplot(model.tseries, [model.numS ,model.numE, model.numI_pre, model.numI_asym, model.numR, model.numQ_S, model.numQ_E, model.numQ_pre, model.numQ_asym, model.numQ_R], labels = ['S', 'E', 'pre', 'I', 'R', 'Q_S', 'Q_E', 'Q_pre', 'Q_I', 'Q_R'], colors = ['darkblue', 'darkred', 'darkgreen', 'indigo', 'darkgoldenrod', 'blue', 'red', 'green', 'blueviolet', 'orange'])
ax2.stackplot(model.tseries, [model.numS ,model.numE, model.numI_pre, model.numI_asym, model.numR, numQ], labels = ['S', 'E', 'pre', 'I', 'R', 'Q'], colors = ['darkblue', 'darkred', 'darkgreen', 'indigo', 'darkgoldenrod', 'brown'])
ax2.vlines([1, 3, 8, 10, 15, 17, 22, 24], 0, numNodes, color = 'black')
my_plot_layout(ax = ax2)
ax2.legend()
plt.show()