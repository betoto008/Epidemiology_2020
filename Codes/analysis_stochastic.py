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

fig, ax = plt.subplots(figsize = (10,8))
p=1.0
sigma=1000

if(p==0.0):
        R0s = np.array([4.5, 3.0, 2.0, 1.2, 0.8])
        Ts = 1- ((meanDegree)/((meanDegree + R0s)))
if(p==1.0):
        R0s = np.array([4.5, 3.0, 2.0, 1.2])

prob_epi_ns = np.empty(100)

#----Load data network of contacts----
N = 2000
graphs_names = np.array(['barabasi-albert','watts-strogatz'])
nodeDegrees = np.loadtxt('../../../../Dropbox/Research/Epidemiology_2020/Text_files/Stochastic/Networks/barabasi-albert/network_degree_distrib_N%d.txt'%(N), dtype=np.int32)
meanDegree = np.mean(nodeDegrees)
meanDegree2 = np.mean(nodeDegrees**2)
degree_distrib = np.histogram(nodeDegrees, bins=range(1, max(nodeDegrees)+1), density = False)
k = degree_distrib[1][:-1]
p_k = degree_distrib[0]/len(nodeDegrees)

ns = np.arange(100)

#----Plot 1----
for i, R0 in enumerate(R0s):
        beta = R0*gamma

        #----Epidemic probability as a function of degree of patient zero----  
        if(p==0.0):
                u_sol = u[np.array([np.sum(p_k*k*(1+(i-1)*T)**(k-1)) for i in u])>(np.sum(p_k*k)*u)][-1]
                print('u:', u_sol)
                prob_epi_n = 1-(np.sum(k*p_k*(1-T+(T*u_sol))**(k-1))/(np.sum(k*p_k)))**ns

        
    prob_epi_n_temp = np.ones(50)
    data = np.loadtxt('../../../../Dropbox/Research/Epidemiology_2020/Text_files/Stochastic/Networks/barabasi-albert/ensemble_I_R0%.1f_sigma%.1f_N%d_p%.1f_barabasi-albert.txt'%(beta/gamma, sigma, N, p))
    max_values = np.array([np.max(data[i,:]) for i in np.arange(len(data[:,0]))])
    data_ext = [((data[i,-1]==0) & (max_values[i] < 20)) for i in np.arange(len(data[:,0]))] #selecting extinct trajectories
    n_ext = len(data[:,0][data_ext])
    prob_ext = n_ext/len(data[:,0])
    data_hist = numpy.histogram(max_values, bins = np.arange(0,int(np.max(max_values))+2, 1), density = True);
    data_hist_ext = numpy.histogram(max_values[data_ext], bins = np.arange(0,int(np.max(max_values))+2, 1), density = True);

    prob_epi_n = 1-(((1-np.cumsum(data_hist_ext[0]))*prob_ext)/(1-np.cumsum(data_hist[0])[data_hist_ext[1][:-1]]))
    prob_epi_ns = np.vstack((prob_epi_ns, prob_epi_n[:80]))
    ax.plot(data_hist_ext[1][:-1]+2, prob_epi_n , 'o', color = colors_R[i], ms = 10,  label = r'$R_0=$%.1f'%(R0))
    #ax.plot(np.array(range(20)), 1-(prob_ext)**(np.array(range(20))),linewidth = 3, color = colors_R[i])
    if(p==1.0):
        ax.plot(np.array(np.arange(20)), 1-(1/R0)**(np.array(np.arange(20))),linewidth = 3, linestyle = '--', color = colors_R[i])
    if(p==0.0):
        ax.plot(np.array(np.arange(20)), e_n,linewidth = 3, linestyle = '--', color = colors_R[i])
ax.hlines(1,0,40)
ax.set_xticks(np.array(np.arange(0,20,2)))
ax.set_xlabel(r'$n$', fontsize = 25)
ax.set_ylabel(r'$P(\mathrm{epi}|n_{\mathrm{max}}\geq n-1)$', fontsize = 25)
ax.set_xlim(0, 15)
ax.set_ylim(-0.05, 1.05)
#ax.set_yscale('log')
ax.tick_params(labelsize = 25)
handles, labels = ax.get_legend_handles_labels()
ax.legend(np.concatenate(([],handles)), np.concatenate(([],labels)) , fontsize = 20, loc = 5, framealpha=.95)
fig.savefig('../Figures/Stochastic/Networks/barabasi-albert/epi_prob_n_SIR_p%.1f.pdf'%(p))

