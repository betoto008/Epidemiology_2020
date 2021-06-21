from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import sys
from models import *
from Epi_models import*
from functions import *
import networkx as nx
import matplotlib.animation as animation
import seaborn as sns
from matplotlib.lines import Line2D
import matplotlib.ticker as ticker
import pickle

Text_files_path = '../../../../Dropbox/Research/Epidemiology_2020/Text_files/'


#----Load data network of contacts----
N = 2000
graphs_names = np.array(['barabasi-albert','watts-strogatz'])
infile_k = open(Text_files_path+'Stochastic/Networks/barabasi-albert/k.pck','rb')
k = pickle.load(infile_k)
infile_k.close()
infile_p_k = open(Text_files_path+'Stochastic/Networks/barabasi-albert/p_k.pck','rb')
p_k = pickle.load(infile_p_k)
infile_p_k.close()

meanDegree = np.sum(k*p_k)
meanDegree2 = np.sum(k**2*p_k)
T_c = meanDegree/(meanDegree2-meanDegree)


ps=[0.0, 1.0]
gamma = 1/6
sigmas=[1000, 1/4]
models = ['SIR', 'SEIR']
R0sp1 = np.flip(np.array([1.2, 2.0, 3.0, 4.5]))
R0sp0 = np.flip(np.array([0.8, 1.2, 2.0, 3.0, 4.5]))
R0S = np.array([R0sp0,R0sp1], dtype="object")



for q, p in enumerate(ps):
	R0s = R0S[q]
	betas= R0s*gamma
	u = np.linspace(0.00005,0.9,100000)
	colors_R = plt.cm.Paired(np.arange(5))

	#----Critical occupancy----
	T_c = meanDegree/(meanDegree2-meanDegree)

	for s, sigma in enumerate(sigmas):
		model = models[s]
		if(sigma==1000):
			Ts = 1-((meanDegree)/((meanDegree+(R0S[q]))))
			u_sols = np.array([u[np.array([np.sum(p_k*k*(1+(i-1)*j)**(k-1)) for i in u])>(np.sum(p_k*k)*u)][-1] for j in Ts])
			S = 1 - np.array([np.sum(p_k*(1+(i-1)*j)**k) for (i, j) in zip(u_sols, Ts)])
		if(sigma==1/4):
			Ts = 1 - ((meanDegree)/((meanDegree + (R0S[q]*gamma)*2*(sigma+gamma)**(-1))))
			u_sols = np.array([u[np.array([np.sum(p_k*k*(1+(i-1)*j)**(k-1)) for i in u])>(np.sum(p_k*k)*u)][-1] for j in Ts])
			S = 1 - np.array([np.sum(p_k*(1+(i-1)*j)**k) for (i, j) in zip(u_sols, Ts)])

		fig2, ax2 = plt.subplots(figsize = (10,8), gridspec_kw={'left':.2})
		fig, ax = plt.subplots(figsize = (10,8), gridspec_kw={'left':.2})
		for b, beta in enumerate(betas):
			u_sol = u_sols[b]
			e_k = 1-(1-Ts[b]+(Ts[b]*u_sol))**k
			Log_Likelihood = np.log10(e_k/S[b])
			ax.plot(k/meanDegree, Log_Likelihood, color = colors_R[b], linestyle = '--', linewidth = 3)


			#----Load data with simulation outcomes----
			data_stats = np.loadtxt('../../../../Dropbox/Research/Epidemiology_2020/Text_files/Stochastic/Networks/barabasi-albert/stats_R0%.1f_sigma%.1f_N%d_p%.1f_barabasi-albert.txt'%(beta/gamma, sigma, N, p))
			data_I = np.loadtxt('../../../../Dropbox/Research/Epidemiology_2020/Text_files/Stochastic/Networks/barabasi-albert/ensemble_I_R0%.1f_sigma%.1f_N%d_p%.1f_barabasi-albert.txt'%(beta/gamma, sigma, N, p))

			max_values = np.array([np.max(data_I[i,:]) for i in np.arange(len(data_I[:,0]))])

			data_ext = np.array([((data_I[i,-1]==0) & (max_values[i] < 20)) for i in np.arange(len(data_I[:,0]))])

			succ_nodes = data_stats[:,0][~data_ext]
			ext_nodes = data_stats[:,0][data_ext]

			z_succ_nodes = len(succ_nodes)
			z_ext_nodes = len(ext_nodes)

			total = len(succ_nodes)+len(ext_nodes)

			### Calculate histograms
			#data = np.histogram(data_stats[:,0], bins=np.logspace(np.log10(2), np.log10(np.max(data_stats[:,0])), 12), density = False)
			#data_succ = np.histogram(np.concatenate((succ_nodes,[])), bins=np.logspace(np.log10(2), np.log10(np.max(data_stats[:,0])), 12), density = False)
			#data_succ_cum = np.cumsum(data_succ[0])/z_succ_nodes
			#data_cum = np.cumsum(data[0])/total
			#Log_Likelihood_data = np.log10(np.gradient(data_succ_cum)/np.gradient(data_cum))[~np.isnan(np.log10(np.gradient(data_succ_cum)/np.gradient(data_cum)))]
			#degrees_data = data[1][:-1][~np.isnan(np.log10(np.gradient(data_succ_cum)/np.gradient(data_cum)))]
			#degrees_data = degrees_data[~np.isinf(Log_Likelihood_data)]
			#Log_Likelihood_data = Log_Likelihood_data[~np.isinf(Log_Likelihood_data)]

			#----Load data with simulation outcomes----
			data_stats = np.loadtxt('../../../../Dropbox/Research/Epidemiology_2020/Text_files/Stochastic/Networks/barabasi-albert/stats_R0%.1f_sigma%.1f_N%d_p%.1f_barabasi-albert.txt'%(beta/gamma, sigma, N, p))
			data_I = np.loadtxt('../../../../Dropbox/Research/Epidemiology_2020/Text_files/Stochastic/Networks/barabasi-albert/ensemble_I_R0%.1f_sigma%.1f_N%d_p%.1f_barabasi-albert.txt'%(beta/gamma, sigma, N, p))
			max_values = np.array([np.max(data_I[i,:]) for i in np.arange(len(data_I[:,0]))])

			#data_ext = np.array([((data_I[i,-1]==0) & (max(data_I[i,:]) < 15)) for i in range(len(data_I[:,0]))])
			data_ext = np.array([((data_I[i,-1]==0) and (max_values[i]<20)) for i in np.arange(len(data_I[:,0]))]) #selecting extinct trajectories
			#data_ext = np.array([((np.max(data_I[i,:])<.5*n_est)) for i in range(len(data_I[:,0]))])
			n_ext = len(data_stats[:,0][data_ext])

			degrees_epi, counts_epi = np.unique(data_stats[:,0][~data_ext], return_counts=True)
			degrees_ext, counts_ext = np.unique(data_stats[:,0][data_ext], return_counts=True)

			max_degree = np.max(np.concatenate((degrees_epi, degrees_ext)))
			degrees = np.array(np.arange(1,int(max_degree)))
			prob_epi_k0_data = np.array([])
			n_epi = 0
			n_ext = 0

			for d in degrees[:]:
			        n_epi = max(0, counts_epi[degrees_epi==d])
			        n_ext = max(0, counts_ext[degrees_ext==d])
			        n_total = n_epi + n_ext
			        if(n_total>0):
			            prob_epi_k0_data = np.append(prob_epi_k0_data, (n_epi)/n_total)
			        else:
			            degrees = degrees[~np.isin(degrees, d)]

			Log_Likelihood_data = np.log10(prob_epi_k0_data/S[b])

			#degrees_array = np.linspace(np.min(degrees_data), np.max(degrees_data), 200)            
			ax.plot(degrees/meanDegree, Log_Likelihood_data,linestyle = '', marker = '.', color=colors_R[b], ms = 15, label = r'$R_0=$%.1f'%(Ts[b]/T_c))

			#ax.plot(degrees_data/meanDegree, Log_Likelihood_data,linestyle = '', marker = '.', color=colors_R[b], ms = 15, label = r'$R_0=$%.1f'%(Ts[b]/T_c))
			
			#ax.plot(k, (slope*meanDegree)*np.log10(k)-(slope*meanDegree)*np.log10(3), color = colors_R[b], linestyle = '--', linewidth = 2)
		ax.hlines(0, 0, 80, linestyle = '--', color = 'silver')
		my_plot_layout(ax = ax, yscale='linear', xscale='log', ylabel=r'$\log{\left(\frac{p(k|\mathrm{epi})}{p(k)}\right)}$', xlabel=r'$k_0/\left\langle k \right\rangle$')
		my_plot_layout(ax = ax2, yscale='linear', xscale='log', xlabel=r'$k_0/\left\langle k \right\rangle$')

		ax.legend(fontsize = 22)
		ax.set_ylim(-.1, .55)
		ax.set_xlim(1.9/meanDegree, 50/meanDegree)
		fig.savefig('../Figures/Stochastic/Networks/barabasi-albert/log-likelihood_'+model+'_p%.1f.pdf'%( p))
		plt.close(fig2)
	    

