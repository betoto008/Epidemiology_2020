#---------Import---------
import sys
import time
from tqdm import tqdm
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

N_ensemble = 200
numNodes = 20
baseGraph    = networkx.barabasi_albert_graph(n=numNodes, m=2)
G_normal     = custom_exponential_graph(baseGraph, scale=100)
# Social distancing interactions:
#G_distancing = custom_exponential_graph(baseGraph, scale=10)
# Quarantine interactions:
#G_quarantine = custom_exponential_graph(baseGraph, scale=5)

fig, ax = plt.subplots(figsize = (10,8))

Q_days = np.array([0, 2, 5, 14])

colors_d = plt.cm.Blues(np.linspace(0,1,len(Q_days)+1))

for i, d in enumerate(Q_days):
	I_cum = np.array([])
	print('RUNNING ENSEMBLE with d = %d days... \n'%(d))
	for n_ensemble in tqdm(np.arange(N_ensemble)):
		model = ExtSEIRSNetworkModel(G=G_normal, beta=0.4, sigma=10000, lamda = 1/4, gamma = 1/6, gamma_asym=1/6, a = 1,  p=1,
		                          G_Q=G_normal, initI_asym=1, isolation_time = d, sigma_Q = 1000, lamda_Q = 1/4, gamma_Q_asym = 1/6,
		                          transition_mode = 'time_in_state', prevalence_ext = 1e-2, eta = 1, mu_H = 1, xi = 1e-10, nu = 1e-10)

		run_tti_sim(model = model, T=7*4, testing_cadence = 'semiweekly2', testing_compliance_random = np.array([True]*model.numNodes), isolation_compliance_positive_individual = [True]*model.numNodes,
			isolation_compliance_positive_groupmate = np.array([True]*model.numNodes), isolation_groups = np.array([range(numNodes)]*model.numNodes), 
			isolation_lag_positive = 0, print_report = False)

		I_cum = np.append(I_cum, (model.numR[-1] + model.numQ_R[-1]))

	ax.hist(I_cum, bins = np.arange(21), density = True, color = colors_d[i+1], label = '%d'%(d))
	ax.vlines(np.mean(I_cum), 0, ax.get_ylim()[1], color = colors_d[i+1])

ax.legend(title = 'Isolations days')

fig.savefig('../../../Figures/Lolli/1_Extended_Model/histograms.pdf')



