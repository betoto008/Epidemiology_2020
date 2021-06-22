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


colors_R = plt.cm.Paired(range(7))
models = ['SEIR', 'SIR']
gamma = 1/6
ps=[0.0]
sigmas=[1/4, 1000]
u = np.linspace(0.00005,0.9,100000)



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

k0s = np.arange(1, 60, 2)

for q, p in enumerate(ps):
        if(p==0.0):
                R0s = np.array([4.5, 3.0, 2.0, 1.2, 0.8])
        if(p==1.0):
                R0s = np.array([4.5, 3.0, 2.0, 1.2])
        for s, sigma in enumerate(sigmas):
                model = models[s]
                #----Plot 1----
                fig, ax = plt.subplots(figsize = (10,8))
                fig2, ax2 = plt.subplots(figsize = (10,8))

                R0_array = np.linspace(1.01, 5, 30)
                T_array = R0_array*T_c
                u_sol_array = np.array([u[np.array([np.sum(p_k*k*(1+(i-1)*T)**(k-1)) for i in u])>(np.sum(p_k*k)*u)][-1] for T in T_array])

                if(sigma==1000):
                        Ts = 1- ((meanDegree)/((meanDegree + R0s)))
                        u_sols = np.array([u[np.array([np.sum(p_k*k*(1+(i-1)*T)**(k-1)) for i in u])>(np.sum(p_k*k)*u)][-1] for T in Ts])
                        if(p==1.0):
                                x, y = np.meshgrid(R0_array, k0s)
                                z = 1-(1/x) 
                        if(p==0.0):
                                x, y = np.meshgrid(R0_array, k0s)
                                x2, y2 = np.meshgrid(u_sol_array, k0s)
                                z = 1-(1-(x*T_c)+((x*T_c)*x2))**y

                if(sigma==1/4):
                        #Ts = 1- ((meanDegree)/((meanDegree + (R0s*gamma)*2*(sigma+gamma)**(-1))))
                        Ts = 1- ((meanDegree)/(meanDegree + np.sqrt(1-4*((sigma*gamma-sigma*((R0s*gamma)))/(sigma+gamma)**2))))
                        u_sols = np.array([u[np.array([np.sum(p_k*k*(1+(i-1)*T)**(k-1)) for i in u])>(np.sum(p_k*k)*u)][-1] for T in Ts])
                        if(p==1.0):
                                x, y = np.meshgrid(R0_array, k0s)
                                z = 1-(1/x)**(2)
                        if(p==0.0):
                                x, y = np.meshgrid(R0_array, k0s)
                                x2, y2 = np.meshgrid(u_sol_array, k0s)
                                z = 1-(1-(x*T_c)+((x*T_c)*x2))**(2*y)

                cs = ax2.contourf(x, y/meanDegree, z, levels = np.linspace(0,1,80), cmap = plt.cm.jet)
                cs2 = ax2.contour(cs, levels=[0.5], colors='k', linestyles = 'dashed', linewidths = 4)

                for r, R0 in enumerate(R0s):

                        beta = R0*gamma

                        #----Epidemic probability as a function of degree of patient zero---- 
                        if(sigma==1000):
                                if(p==0.0):
                                        R0 = Ts[r]/T_c
                                        u_sol = u_sols[r]
                                        prob_epi_k0 = 1-(1-Ts[r]+(Ts[r]*u_sol))**k0s
                                if(p==1.0):
                                        prob_epi_k0 = np.ones_like(k0s)*(1-(1/R0))

                        if(sigma==1/4):
                                if(p==0.0):
                                        R0 = Ts[r]/T_c
                                        u_sol = u_sols[r]
                                        prob_epi_k0 = 1-(1-Ts[r]+(Ts[r]*u_sol))**(2*k0s)
                                if(p==1.0):
                                        R0 = np.sqrt(1-4*((sigma*gamma-sigma*beta)/(sigma+gamma)**2))
                                        prob_epi_k0 = np.ones_like(k0s)*(1-(1/R0)**2)
                        

                        #----Load data with simulation outcomes----
                        data_nodes = np.loadtxt('../../../../Dropbox/Research/Epidemiology_2020/Text_files/Stochastic/Networks/barabasi-albert/stats_R0%.1f_sigma%.1f_N%d_p%.1f_barabasi-albert.txt'%(beta/gamma, sigma, N, p))
                        data_I = np.loadtxt('../../../../Dropbox/Research/Epidemiology_2020/Text_files/Stochastic/Networks/barabasi-albert/ensemble_I_R0%.1f_sigma%.1f_N%d_p%.1f_barabasi-albert.txt'%(beta/gamma, sigma, N, p))
                        data_E = np.loadtxt('../../../../Dropbox/Research/Epidemiology_2020/Text_files/Stochastic/Networks/barabasi-albert/ensemble_E_R0%.1f_sigma%.1f_N%d_p%.1f_barabasi-albert.txt'%(beta/gamma, sigma, N, p))
                        nodes_ext , nodes_epi, I_ext, I_succ, max_values, data_ext = sort_nodes(p=p, beta=beta, sigma=sigma, gamma=gamma, data_I = data_I, data_E = data_E, data_nodes = data_nodes[:,0], upper_limit = 2)

                        degrees_epi, counts_epi = np.unique(nodes_epi, return_counts=True)
                        degrees_ext, counts_ext = np.unique(nodes_ext, return_counts=True)

                        max_degree = np.max(np.concatenate((degrees_epi, degrees_ext)))
                        degrees = np.array(np.arange(1,int(max_degree)))
                        prob_epi_k0_data = np.array([])
                        n_epi = 0
                        n_ext = 0

                        for d in degrees[:]:
                                if(counts_ext[degrees_ext==d].size>0 and counts_epi[degrees_epi==d].size>0):
                                        n_epi = counts_epi[degrees_epi==d][0]
                                        #print(n_epi)
                                        n_ext = counts_ext[degrees_ext==d][0]
                                        n_total = n_epi + n_ext
                                        prob_epi_k0_data = np.append(prob_epi_k0_data, (n_epi)/n_total)
                                        print(d, n_total)
                                else:
                                        degrees = degrees[~np.isin(degrees, d)]
                                
                        
                        ax.plot(degrees, prob_epi_k0_data , '^', color = colors_R[r], ms = 12,  label = r'$R_0=$%.1f'%(R0))
                        ax.plot(k0s, prob_epi_k0,linewidth = 2, linestyle = '-', color = colors_R[r])

                        for j in np.array([int(i) for i in np.logspace(0, np.log10(len(prob_epi_k0_data)-1), 8)]):
                                ax2.scatter(R0, (degrees[j])/meanDegree, marker = 's', color = plt.cm.jet(np.linspace(0,1,80))[int(79*prob_epi_k0_data[j])], s = 200, edgecolors='k')

                # Plot 1
                ax.hlines(1,1,40, linestyle = '--', color = 'silver')
                my_plot_layout(ax = ax, xlabel = r'$k_0/\left\langle k \right\rangle$', ylabel = r'$P(\mathrm{epi}|k_0)$', yscale = 'linear', xscale='linear')
                ax.set_xticks(np.array(np.arange(1,40,5)))
                ax.set_xlim(1, 40)
                ax.set_ylim(-0.05, 1.05)
                handles, labels = ax.get_legend_handles_labels()
                ax.legend(np.concatenate(([],handles)), np.concatenate(([],labels)) , fontsize = 20, loc = 4, framealpha=.95)
                fig.savefig('../Figures/Stochastic/Networks/barabasi-albert/Epi_prob_k0_'+model+'_p%.1f.png'%(p))

                # Plot 2
                my_plot_layout(ax=ax2, xlabel=r'$R_0$', ylabel=r'$k_0/\left\langle k \right\rangle$', yscale='log')
                if(p==1.0):
                        if(sigma==1000):
                                ax2.set_xlim(1.02, 4.6)
                        if(sigma==1/4):
                                ax2.set_xlim(1.02, 2.3)
                if(p==0.0):
                        ax2.set_xlim(1.02, 4)
                ax2.set_ylim(1.9/meanDegree, 50/meanDegree)
                cbar = fig2.colorbar(cs, ticks=np.linspace(0,1,5))
                #cbar.set_ticks(np.arange(5))
                cbar.set_label('Probability of epidemic', fontsize = 25)
                cbar.set_ticklabels(np.linspace(0,1,5))
                cbar.ax.tick_params(labelsize = 25)
                cbar.add_lines(cs2)
                fig2.savefig('../Figures/Stochastic/Networks/barabasi-albert/Epi_prob_k0_'+model+'_p%.1f_HM.png'%(p))
                


