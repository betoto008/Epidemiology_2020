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
from matplotlib.ticker import FormatStrFormatter
import matplotlib.ticker as mtick

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

sample_sizes = np.array([int((N/100)*1.875), int((N/100)*3.125), int((N/100)*5.0)])
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
		fig, ax = plt.subplots(figsize=(12,8))
		fig2, ax2 = plt.subplots(figsize=(12,8))
		fig3, ax3 = plt.subplots(figsize=(12,8))
		fig4, ax4 = plt.subplots(figsize=(12,8))
		fig5, ax5 = plt.subplots(figsize=(12,8))
		fig6, ax6 = plt.subplots(figsize=(12,8))

		avg_t_uniform = np.array([])
		avg_cs_uniform = np.array([])
		avg_t_aposteriori = np.array([])
		avg_cs_aposteriori = np.array([])

		for beta, T_total in zip(betas, Ts_total):
	        
			print('p:', p, 'sig:', sigma, 'R0:', beta/gamma, 'T:', T_total)
	
			for m, color, color_m in zip(sample_sizes, colors, colors_m):
				data_sampling_uniform = np.loadtxt(Text_files_path+'Sampling/Networks/barabasi-albert/uniform/sampling_stats_R0%.1f_sigma%.1f_N%d_p%.1f_m%d_barabasi-albert.txt'%(beta/gamma,sigma,N,p,m))
				data_sampling_uniform2 = np.loadtxt(Text_files_path+'Sampling/Networks/barabasi-albert/uniform/prueba_CS/sampling_stats_R0%.1f_sigma%.1f_N%d_p%.1f_m%d_barabasi-albert.txt'%(beta/gamma,sigma,N,p,m))
				data_sampling_uniform_times = np.concatenate((data_sampling_uniform[data_sampling_uniform[:,0]!=1000,0],data_sampling_uniform2[data_sampling_uniform2[:,0]!=1000,0]))
				data_sampling_uniform_cs = data_sampling_uniform2[data_sampling_uniform2[:,0]!=1000,5]

				data_sampling_aposteriori = np.loadtxt(Text_files_path+'Sampling/Networks/barabasi-albert/aposteriori/sampling_stats_R0%.1f_sigma%.1f_N%d_p%.1f_m%d_barabasi-albert.txt'%(beta/gamma,sigma,N,p,m))
				data_sampling_aposteriori2 = np.loadtxt(Text_files_path+'Sampling/Networks/barabasi-albert/aposteriori/prueba_CS/sampling_stats_R0%.1f_sigma%.1f_N%d_p%.1f_m%d_barabasi-albert.txt'%(beta/gamma,sigma,N,p,m))
				data_sampling_aposteriori_times = np.concatenate((data_sampling_aposteriori[data_sampling_aposteriori[:,0]!=1000,0],data_sampling_aposteriori2[data_sampling_aposteriori2[:,0]!=1000,0]))
				data_sampling_aposteriori_cs = data_sampling_aposteriori2[data_sampling_aposteriori2[:,0]!=1000,5]


				avg_t_uniform = np.append(avg_t_uniform, np.mean(data_sampling_uniform_times))
				avg_cs_uniform = np.append(avg_cs_uniform, np.mean(data_sampling_uniform_cs))

				avg_t_aposteriori = np.append(avg_t_aposteriori, np.mean(data_sampling_aposteriori_times))
				avg_cs_aposteriori = np.append(avg_cs_aposteriori, np.mean(data_sampling_aposteriori_cs))
			i_b += 1

		avg_t_uniform = np.reshape(avg_t_uniform, (np.size(sample_sizes), np.size(Ts_total)))
		avg_cs_uniform = np.reshape(avg_cs_uniform, (np.size(sample_sizes), np.size(Ts_total)))
		avg_t_aposteriori = np.reshape(avg_t_aposteriori, (np.size(sample_sizes), np.size(Ts_total)))
		avg_cs_aposteriori = np.reshape(avg_cs_aposteriori, (np.size(sample_sizes), np.size(Ts_total)))

		sns.heatmap(avg_t_uniform, ax = ax, cmap=plt.cm.twilight, center = 0, cbar = True, cbar_kws={'label': r'$T_U$'})
		sns.heatmap(avg_cs_uniform/meanDegree, ax = ax2, cmap=plt.cm.seismic, center = 1, cbar = True, cbar_kws={'label': r'$n_U/\left\langle k \right\rangle$'})
		sns.heatmap(avg_t_aposteriori, ax = ax3, cmap=plt.cm.twilight, center = 0, cbar = True, cbar_kws={'label': r'$T_A$'})
		sns.heatmap(avg_cs_aposteriori/meanDegree, ax = ax4, cmap=plt.cm.seismic, center = 1, cbar = True, cbar_kws={'label': r'$n_A/\left\langle k \right\rangle$'})
		sns.heatmap(avg_t_aposteriori/avg_t_uniform, ax = ax5, cmap=plt.cm.seismic, center = 1, cbar = True, cbar_kws={'label': r'$T_{A}/T_{U}$'})
		sns.heatmap(avg_cs_aposteriori/avg_cs_uniform, ax = ax6, cmap=plt.cm.seismic, center = 1, cbar = True, cbar_kws={'label': r'$n_{A}/n_{U}$'})



		my_plot_layout(ax=ax, ylabel=r'Sample size (%)', xlabel=r'$R_0$', yscale='linear')
		ax.set_yticks([.5, 1.5, 2.5])
		ax3.set_yticklabels(FormatStrFormatter('%.2f').format_ticks(sample_sizes/N*100))
		ax.set_xticks(np.array([g + 0.5 for g in np.arange(len(betas))]))
		ax.set_xticklabels(betas/gamma)
		cbar = ax.collections[0].colorbar
		cbar.ax.tick_params(labelsize=18)
		ax.figure.axes[-1].yaxis.label.set_size(24)
		fig.savefig('../Figures/Sampling/Networks/barabasi-albert/prueba/avg_t_uniform_'+model+'_p%.1f.png'%(p))
		plt.close(fig)

		my_plot_layout(ax=ax2, ylabel=r'Sample size (%)', xlabel=r'$R_0$', yscale='linear')
		ax2.set_yticks([.5, 1.5, 2.5])
		ax2.set_yticklabels(FormatStrFormatter('%.2f').format_ticks(sample_sizes/N*100))
		ax2.set_xticks(np.array([g + 0.5 for g in np.arange(len(betas))]))
		ax2.set_xticklabels(betas/gamma)
		cbar = ax2.collections[0].colorbar
		cbar.ax.tick_params(labelsize=18)
		ax2.figure.axes[-1].yaxis.label.set_size(24)
		fig2.savefig('../Figures/Sampling/Networks/barabasi-albert/prueba/avg_cs_uniform_'+model+'_p%.1f.png'%(p))
		plt.close(fig2)

		my_plot_layout(ax=ax3, ylabel=r'Sample size (%)', xlabel=r'$R_0$', yscale='linear')
		ax3.set_yticks([.5, 1.5, 2.5])
		ax3.set_yticklabels(FormatStrFormatter('%.2f').format_ticks(sample_sizes/N*100))
		ax3.set_xticks(np.array([g + 0.5 for g in np.arange(len(betas))]))
		ax3.set_xticklabels(betas/gamma)
		cbar = ax3.collections[0].colorbar
		cbar.ax.tick_params(labelsize=18)
		ax3.figure.axes[-1].yaxis.label.set_size(24)
		fig3.savefig('../Figures/Sampling/Networks/barabasi-albert/prueba/avg_t_aposteriori_'+model+'_p%.1f.png'%(p))
		plt.close(fig3)

		my_plot_layout(ax=ax4, ylabel=r'Sample size (%)', xlabel=r'$R_0$', yscale='linear')
		ax4.set_yticks([.5, 1.5, 2.5])
		ax4.set_yticklabels(FormatStrFormatter('%.2f').format_ticks(sample_sizes/N*100))
		ax4.set_xticks(np.array([g + 0.5 for g in np.arange(len(betas))]))
		ax4.set_xticklabels(betas/gamma)
		cbar = ax4.collections[0].colorbar
		cbar.ax.tick_params(labelsize=18)
		ax4.figure.axes[-1].yaxis.label.set_size(24)
		fig4.savefig('../Figures/Sampling/Networks/barabasi-albert/prueba/avg_cs_aposteriori_'+model+'_p%.1f.png'%(p))
		plt.close(fig4)

		my_plot_layout(ax=ax5, ylabel=r'Sample size (%)', xlabel=r'$R_0$', yscale='linear')
		ax5.set_yticks([.5, 1.5, 2.5])
		ax5.set_yticklabels(FormatStrFormatter('%.2f').format_ticks(sample_sizes/N*100))
		ax5.set_xticks(np.array([g + 0.5 for g in np.arange(len(betas))]))
		ax5.set_xticklabels(betas/gamma)
		cbar = ax5.collections[0].colorbar
		cbar.ax.tick_params(labelsize=18)
		ax5.figure.axes[-1].yaxis.label.set_size(24)
		fig5.savefig('../Figures/Sampling/Networks/barabasi-albert/prueba/figure_times_'+model+'_p%.1f.png'%(p))
		plt.close(fig5)

		my_plot_layout(ax=ax6, ylabel=r'Sample size (%)', xlabel=r'$R_0$', yscale='linear')
		ax6.set_yticks([.5, 1.5, 2.5])
		ax6.set_yticklabels(FormatStrFormatter('%.2f').format_ticks(sample_sizes/N*100))
		ax6.set_xticks(np.array([g + 0.5 for g in np.arange(len(betas))]))
		ax6.set_xticklabels(betas/gamma)
		cbar = ax6.collections[0].colorbar
		cbar.ax.tick_params(labelsize=18)
		ax6.figure.axes[-1].yaxis.label.set_size(24)
		fig6.savefig('../Figures/Sampling/Networks/barabasi-albert/prueba/figure_cs_'+model+'_p%.1f.png'%(p))
		plt.close(fig6)

		i_p+=1

	        
	        
	        


