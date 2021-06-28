from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import sys
from models import *
from Epi_models import*
from functions import *
import networkx as nx
import matplotlib.animation as animation
import seaborn
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

sample_sizes = [int((N/100)*1.875), int((N/100)*3.125), int((N/100)*5.0)]
gamma = 1/6
sigmas=[1/4, 1000]
models = ['SEIR', 'SIR']
ps = np.array([0.0, 1.0])
R0s = np.array([1.2, 2.0, 3.0, 4.5])
R0s2 = np.array([0.8, 1.2, 2.0, 3.0, 4.5])
colors_m = ['darkgreen', 'darkblue', 'darkred']
colors = ['green', 'blue', 'red']
i_p = 0

for s, sigma in enumerate(sigmas):

	Betas = np.array([R0s2*gamma,R0s*gamma], dtype=object)
	#print('Betas:', Betas)
	R0_Ns = (R0s2*(meanDegree2-meanDegree))/(meanDegree*(meanDegree+R0s2))
	R0_Es = np.sqrt(1-4*((sigma*gamma-sigma*Betas[0])/(sigma + gamma)**2))
	R0_ENs = (R0_Es*(meanDegree2-meanDegree))/(meanDegree*(meanDegree+R0_Es))

	#print('R0s:', R0s2)
	#print('R0_Es:',R0_Es)
	#print('R0_Ns:',R0_Ns)
	#print('R0_ENs:',R0_ENs)

	lambdas = ((-sigma-gamma)/(2)) + (1/2)*np.sqrt((sigma-gamma)**2 + 4*sigma*Betas[1]) #exponential growth rates

	ests = (1/(lambdas)).astype(int)
	if(sigma==1/4):
	    est_Ns = (1/((R0_ENs-1)*gamma)).astype(int)
	if(sigma==1000):
	    est_Ns = (1/((R0_Ns-1)*gamma)).astype(int)

	Ts_Total = np.array([5*est_Ns,5*ests], dtype="object") #for R0 = [1.2, 1.5, 3.0, 4.5] use [160, 35, 20] and half of it if when sigma=1000
	print('Times:',Ts_Total)

	model = models[s]
	for p, Ts_total, betas in zip(ps, Ts_Total, Betas):
		i_b = 0
		for beta, T_total in zip(betas, Ts_total):
			fig, ax = plt.subplots(figsize=(12,8))
			fig2, ax2 = plt.subplots(figsize=(12,8))
			print('p:', p, 'sig:', sigma, 'R0:', beta/gamma, 'T:', T_total)
			days = np.arange(T_total+1)
			#data_I = np.loadtxt(Text_files_path+'Stochastic/Networks/barabasi-albert/ensemble_I_R0%.1f_sigma%.1f_N%d_p%.1f_barabasi-albert.txt'%(beta/gamma,sigma,N,p))
			#I_avg = np.mean(data_I , axis = 0)
			#max_values = np.array([np.max(data_I[i,:]) for i in np.arange(len(data_I[:,0]))])
			#data_ext = np.array([((data_I[i,-1]==0) and (max_values[i]<20)) for i in np.arange(len(data_I[:,0]))]) #selecting extinct trajectories
			#data_I2 = data_I[~data_ext,:]
			#I_avg_not_ext = np.mean(data_I2 , axis = 0)
			#I_det = np.loadtxt(Text_files_path+'Deterministic/deterministic_R0%.1f_sigma%.1f_p%.1f_N%d.txt'%(beta/gamma, sigma, p, N))

			#ax = plot_cum_prob_time(T_total, sample_sizes , beta/gamma, sigma, p, N, days , I_avg_not_ext, I_avg_not_ext, colors, net_name = graphs_names[0], 
			                                #folder=Text_files_path+'Sampling/Networks/'+graphs_names[0], external_ax = ax)
			ax.hlines([0, .9, 1],0, T_total, linestyle = 'dashed', color = 'silver', linewidth = 3)
			#ax2.hlines([0, .9, 1],0,T_total, linestyle = 'dashed', color = 'silver', linewidth = 3)
			for m, color, color_m in zip(sample_sizes, colors, colors_m):
				data_sampling_uniform = np.loadtxt(Text_files_path+'Sampling/Networks/barabasi-albert/uniform/k_normalization/sampling_stats_R0%.1f_sigma%.1f_N%d_p%.1f_m%d_barabasi-albert.txt'%(beta/gamma,sigma,N,p,m))
				#data_sampling_uniform2 = np.loadtxt(Text_files_path+'Sampling/Networks/barabasi-albert/uniform/prueba_CS/sampling_stats_R0%.1f_sigma%.1f_N%d_p%.1f_m%d_barabasi-albert.txt'%(beta/gamma,sigma,N,p,m))
				#data_sampling_uniform_times = np.concatenate((data_sampling_uniform[data_sampling_uniform[:,0]!=1000,0],data_sampling_uniform2[data_sampling_uniform2[:,0]!=1000,0]))
				data_sampling_uniform_times = data_sampling_uniform[data_sampling_uniform[:,0]!=1000,0]
				data_sampling_uniform_cs = data_sampling_uniform[data_sampling_uniform[:,0]!=1000,5]

				data_sampling_aposteriori = np.loadtxt(Text_files_path+'Sampling/Networks/barabasi-albert/aposteriori/k_normalization/sampling_stats_R0%.1f_sigma%.1f_N%d_p%.1f_m%d_barabasi-albert.txt'%(beta/gamma,sigma,N,p,m))
				#data_sampling_aposteriori2 = np.loadtxt(Text_files_path+'Sampling/Networks/barabasi-albert/aposteriori/prueba_CS/sampling_stats_R0%.1f_sigma%.1f_N%d_p%.1f_m%d_barabasi-albert.txt'%(beta/gamma,sigma,N,p,m))
				#data_sampling_aposteriori_times = np.concatenate((data_sampling_aposteriori[data_sampling_aposteriori[:,0]!=1000,0],data_sampling_aposteriori2[data_sampling_aposteriori2[:,0]!=1000,0]))
				data_sampling_aposteriori_times = data_sampling_aposteriori[data_sampling_aposteriori[:,0]!=1000,0]
				data_sampling_aposteriori_cs = data_sampling_aposteriori[data_sampling_aposteriori[:,0]!=1000,5]


				data_uniform_times = np.histogram(data_sampling_uniform_times, bins = np.arange(T_total+1), density=True)
				data_aposteriori_times = np.histogram(data_sampling_aposteriori_times, bins = np.arange(T_total+1), density=True)
				data_uniform_cs = np.histogram(data_sampling_uniform_cs, bins = np.arange(1, 100), density=True)
				data_aposteriori_cs = np.histogram(data_sampling_aposteriori_cs, bins = np.arange(1, 100), density=True)


				ax.plot(np.arange(T_total), np.cumsum(data_uniform_times[0]), linestyle = '--', marker = 'o', color = color, ms = 6, linewidth = 1, label = '$m=%.1f %%$'%(m*100/N))
				ax.plot(np.arange(T_total), np.cumsum(data_aposteriori_times[0]), linestyle = '--', marker = 's', color = color, ms = 6, linewidth = 1)

				ax2.plot(data_uniform_cs[1][:-1], data_uniform_cs[0], linestyle = '--', marker = 'o', color = color, ms = 6, linewidth = 1, label = '$m=%.1f %%$'%(m*100/N))
				ax2.plot(data_aposteriori_cs[1][:-1], data_aposteriori_cs[0], linestyle = '--', marker = 's', color = color, ms = 6, linewidth = 1)

			my_plot_layout(ax=ax, xlabel='Time [days]', ylabel='Cum. Prob. detection', yscale='linear')
			ax.set_xlim(0, np.min((90, T_total)))
			ax.legend(loc = 4, fontsize = 26)
			fig.savefig('../Figures/Sampling/Networks/barabasi-albert/prueba/prob_cum_time_'+model+'_R0%.1f_p%.1f.png'%(beta/gamma,p))
			plt.close(fig)

			my_plot_layout(ax=ax2, xlabel=r'Cluster size $n$', ylabel='Prob. detection', yscale='linear')
			ax2.set_xlim(0, 100)
			#ax2.set_ylim(0, np.max(data_aposteriori[0]*2.5))
			ax2.legend(loc = 4, fontsize = 26)
			fig2.savefig('../Figures/Sampling/Networks/barabasi-albert/prueba/prob_cs_'+model+'_R0%.1f_p%.1f.png'%(beta/gamma,p))
			plt.close(fig2)

			i_b += 1
	        
	        


