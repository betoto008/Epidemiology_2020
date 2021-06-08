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

Text_files_path = '../../../../Dropbox/Research/Epidemiology_2020/Text_files/'


colors_R = plt.cm.Paired(range(7))
models = ['SIR', 'SEIR']
gamma = 1/6
ps=[1.0, 0.0]
sigmas=[1000, 1/4]
u = np.linspace(0.00005,0.9,100000)



#----Load data network of contacts----
N = 2000
graphs_names = np.array(['barabasi-albert','watts-strogatz'])
nodeDegrees = np.loadtxt('../../../../Dropbox/Research/Epidemiology_2020/Text_files/Stochastic/Networks/barabasi-albert/network_degree_distrib_N%d.txt'%(N), dtype=np.int32)
meanDegree = np.mean(nodeDegrees)
meanDegree2 = np.mean(nodeDegrees**2)
degree_distrib = np.histogram(nodeDegrees, bins=range(1, max(nodeDegrees)+1), density = False)
k = degree_distrib[1][:-1]
p_k = degree_distrib[0]/len(nodeDegrees)
T_c = meanDegree/(meanDegree2-meanDegree)

ns = np.arange(80)

for q, p in enumerate(ps):
        if(p==0.0):
                R0s = np.array([4.5, 3.0, 2.0, 1.2, 0.8])
                Ts = 1- ((meanDegree)/((meanDegree + R0s)))
                u_sols = np.array([u[np.array([np.sum(p_k*k*(1+(i-1)*T)**(k-1)) for i in u])>(np.sum(p_k*k)*u)][-1] for T in Ts])
        if(p==1.0):
                R0s = np.array([4.5, 3.0, 2.0, 1.2])
        for s, sigma in enumerate(sigmas):
                model = models[s]
                #----Plot 1----
                fig, ax = plt.subplots(figsize = (10,8))
                fig2, ax2 = plt.subplots(figsize = (10,8))

                R0_array = np.linspace(1.01, 6, 30)
                T_array = R0_array*T_c
                u_sol_array = np.array([u[np.array([np.sum(p_k*k*(1+(i-1)*T)**(k-1)) for i in u])>(np.sum(p_k*k)*u)][-1] for T in T_array])

                if(sigma==1000):
                        if(p==1.0):
                                x, y = np.meshgrid(R0_array, ns)
                                z = 1-(1/x)**(y) 
                        if(p==0.0):
                                z0 = np.array([np.sum(k*p_k*(1-j+(j*i))**(k-1))/(np.sum(p_k*k)) for (i, j) in zip(u_sol_array, T_array)])
                                x, y = np.meshgrid(R0_array, ns)
                                x2, y2 = np.meshgrid(z0, ns)
                                z = 1-(x2**(y))

                if(sigma==1/4):

                        if(p==1.0):
                                x, y = np.meshgrid(R0_array, ns)
                                z = 1-(1/x)**(2*y)
                        if(p==0.0):
                                z0 = np.array([np.sum(k*p_k*(1-j+(j*i))**(k-1))/(np.sum(p_k*k)) for (i, j) in zip(u_sol_array, T_array)])
                                x, y = np.meshgrid(R0_array, ns)
                                x2, y2 = np.meshgrid(z0, ns)
                                z = 1-(x2**(y))

                cs = ax2.contourf(x, y/meanDegree, z, levels = np.linspace(0,1,20), cmap = plt.cm.jet)
                cs2 = ax2.contour(cs, levels=[0.75], colors='k', linestyles = 'dashed', linewidths = 4)

                for r, R0 in enumerate(R0s):

                        beta = R0*gamma

                        #----Epidemic probability as a function of degree of patient zero---- 
                        if(sigma==1000):
                                if(p==0.0):
                                        R0 = Ts[r]/T_c
                                        u_sol = u_sols[r]
                                        prob_epi_n = 1-(np.sum(k*p_k*(1-Ts[r]+(Ts[r]*u_sol))**(k-1))/(np.sum(k*p_k)))**ns
                                if(p==1.0):
                                        prob_epi_n = 1-(1/R0)**ns

                        if(sigma==1/4):
                                if(p==0.0):
                                        R0 = Ts[r]/T_c
                                        u_sol = u[np.array([np.sum(p_k*k*(1+(i-1)*Ts[r])**(k-1)) for i in u])>(np.sum(p_k*k)*u)][-1]
                                        prob_epi_n = 1-(np.sum(k*p_k*(1-Ts[r]+(Ts[r]*u_sol))**(k-1))/(np.sum(k*p_k)))**(ns)
                                if(p==1.0):
                                        R0 = np.sqrt(1-4*(((1/4)*gamma-sigmas[1]*beta)/((1/4)+gamma)**2))
                                        prob_epi_n = 1-(1/R0)**(2*ns)
 
                        data = np.loadtxt('../../../../Dropbox/Research/Epidemiology_2020/Text_files/Stochastic/Networks/barabasi-albert/ensemble_I_R0%.1f_sigma%.1f_N%d_p%.1f_barabasi-albert.txt'%(beta/gamma, sigma, N, p))
                        max_values = np.array([np.max(data[i,:]) for i in np.arange(len(data[:,0]))])
                        data_ext = np.array([((data[i,-1]==0) & (max_values[i] < 20)) for i in np.arange(len(data[:,0]))]) #selecting extinct trajectories
                        n_ext = len(data[:,0][data_ext])
                        prob_ext = n_ext/len(data[:,0])
                        hist = numpy.histogram(max_values, bins = np.arange(1,int(np.max(max_values))+2), density = True);
                        hist_ext = numpy.histogram(max_values[data_ext], bins = np.arange(1,int(np.max(max_values))+2), density = True);

                        prob_epi_n_data = 1-(((1-np.cumsum(hist_ext[0]))*prob_ext)/(1-np.cumsum(hist[0])))
                        ax.plot((hist_ext[1][:-1]+1)/meanDegree, prob_epi_n_data , '^', color = colors_R[r], ms = 12,  label = r'$R_0=$%.1f'%(R0))
                        ax.plot(ns/meanDegree, prob_epi_n,linewidth = 3, linestyle = '--', color = colors_R[r])

                        for j in np.array([int(i) for i in np.logspace(np.log10(2), np.log10(50), 8)-2]):
                                ax2.scatter(R0, (ns[j]+1)/meanDegree, marker = 's', color = plt.cm.jet(np.linspace(0,1,20))[int(19*prob_epi_n_data[j])], s = 200, edgecolors='k')

                # Plot 1
                ax.hlines(1,0,40)
                ax.set_xticks(np.array(np.arange(0,20,2)))
                ax.set_xlim(0, 15/meanDegree)
                ax.set_ylim(-0.05, 1.05)
                my_plot_layout(ax = ax, xlabel = r'$n/\left\langle k \right\rangle$', ylabel = r'$P(\mathrm{epi}|n_{\mathrm{max}}\geq n-1)$', yscale = 'linear')
                handles, labels = ax.get_legend_handles_labels()
                ax.legend(np.concatenate(([],handles)), np.concatenate(([],labels)) , fontsize = 20, loc = 5, framealpha=.95)
                fig.savefig('../Figures/Stochastic/Networks/barabasi-albert/Epi_prob_n_'+model+'_p%.1f.pdf'%(p))

                # Plot 2
                my_plot_layout(ax=ax2, xlabel=r'$R_0$', ylabel=r'$n/\left\langle k \right\rangle$', yscale='log')
                if(p==1.0):
                        ax2.set_xlim(1.02, 2.3)
                if(p==0.0):
                        ax2.set_xlim(1.02, 4)
                ax2.set_ylim(1.9/meanDegree, 30/meanDegree)
                cbar = fig2.colorbar(cs, ticks=np.linspace(0,1,5))
                cbar.ax.set_ylabel('Probability of epidemic', fontsize = 25)
                cbar.set_ticklabels(np.linspace(0,1,5))
                cbar.ax.tick_params(labelsize = 25)
                cbar.add_lines(cs2)

                fig2.savefig('../Figures/Stochastic/Networks/barabasi-albert/Epi_prob_n_'+model+'_p%.1f_HM.pdf'%(p))
                


