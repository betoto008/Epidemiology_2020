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
models = ['SEIR']
gamma = 1/6
ps=[0.0]
sigmas=[1/4]
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

Ts_array = np.linspace(T_c, T_c*6, 20)
u_sol_array = np.array([u[np.array([np.sum(p_k*k*(1+(i-1)*j)**(k-1)) for i in u])>(np.sum(p_k*k)*u)][-1] for j in Ts_array])
S_SIR = 1 - np.array([np.sum(p_k*(1+(i-1)*j)**k) for (i, j) in zip(u_sol_array, Ts_array)])
S_SEIR = 1 - np.array([np.sum(p_k*(1-j+i*j)**(2*k)) for (i, j) in zip(u_sol_array, Ts_array)])


ps = [1,0]
sigmas = [1000, 1/4]
R0_array = np.linspace(1, 8, 50)
R0s_p1 = np.array([4.5, 3.0, 2.0, 1.2])
R0s_p0 = np.array([4.5, 3.0, 2.0, 1.2, 0.8])
markersR0 = ['*', '^', 'o', 's', 'X']
betas1 = R0s_p1*gamma
betas2 = R0s_p0*gamma
betas_p = np.array([betas1, betas2], dtype = object)
R0s_SIR_p = np.array([R0s_p1, (1-((4)/((4+((betas2)/(gamma))))))/T_c], dtype = object)
#R0s_SEIR_p = np.array([np.sqrt(1-4*(((1/4)*gamma-(1/4)*betas1)/((1/4)+gamma)**2)), (1-((4)/(4+((betas2)*2*(gamma+(1/4))**(-1)))))/T_c])
R0s_SEIR_p = np.array([np.sqrt(1-4*(((1/4)*gamma-(1/4)*betas1)/((1/4)+gamma)**2)), (1-((4)/(4+(np.sqrt(1-4*(((1/4)*gamma-(1/4)*betas2)/((1/4)+gamma)**2))))))/T_c], dtype = object)
R0s_all = np.array([R0s_SIR_p, R0s_SEIR_p], dtype = object)

colors_R = plt.cm.Paired(range(7))
colors_p = ['darkblue', 'darkred']
lines_symbols = []
labels_symbols = ['Simulation']
labels_model = ['SIR', 'SEIR']

for l, sigma in enumerate(sigmas):
    fig, ax = plt.subplots(figsize = (10,8), gridspec_kw={'left':.2})
    if(sigma==1/4):
        ax.set_title('SEIR model', fontsize = 28)
    if(sigma==1000):
        ax.set_title('SIR model', fontsize = 28)
    R0s_p = R0s_all[l]
    for p, betas, R0s, color_p in zip(ps, betas_p, R0s_p, colors_p):
        
        p_epi_array = np.array([])
        i_b = 0
        for beta, R0, color in zip(betas, R0s, colors_R):
            #----Edge Occupancy probability----

            tau = 1/gamma
            data_nodes = np.loadtxt('../../../../Dropbox/Research/Epidemiology_2020/Text_files/Stochastic/Networks/barabasi-albert/stats_R0%.1f_sigma%.1f_N%d_p%.1f_barabasi-albert.txt'%(beta/gamma, sigma, N, p))
            data_I = np.loadtxt('../../../../Dropbox/Research/Epidemiology_2020/Text_files/Stochastic/Networks/barabasi-albert/ensemble_I_R0%.1f_sigma%.1f_N%d_p%.1f_barabasi-albert.txt'%(beta/gamma, sigma, N, p))
            data_E = np.loadtxt('../../../../Dropbox/Research/Epidemiology_2020/Text_files/Stochastic/Networks/barabasi-albert/ensemble_E_R0%.1f_sigma%.1f_N%d_p%.1f_barabasi-albert.txt'%(beta/gamma, sigma, N, p))
            nodes_ext , nodes_succ, I_ext, I_succ, max_values, data_ext = sort_nodes(p=p, beta=beta, sigma=sigma, gamma=gamma, data_I = data_I, data_E = data_E, data_nodes = data_nodes[:,0], upper_limit = 3)
            
            if(sigma==1/4):
                tau = 2*((gamma)+(sigma))**(-1)
            s = ((R0)-1)
            
            #n_epi = len(max_values[max_values>(1/(s/tau))])
            #n_ext = len(max_values[max_values<(1/(s/tau))])
            
            n_epi = len(nodes_succ)
            n_ext = len(nodes_ext)
            
            n_total = n_epi + n_ext
            p_epi_array = np.append(p_epi_array, (n_epi)/n_total)
            
            ax.plot(R0, p_epi_array[-1], marker = markersR0[i_b], color = color_p, ms = 15, linestyle = '')
            i_b +=1
        
        if(p==1.0):
            if(sigma==1000):
                ax.plot(R0_array, 1-((1/R0_array)),linestyle = 'dashed', linewidth = 3, color = color_p, label = 'Fully-connected')
            if(sigma==1/4):
                ax.plot(R0_array, 1-((1/R0_array)**2),linestyle = 'dashed', linewidth = 3, color = color_p, label = 'Fully-connected')
        if(p==0.0):
            if(sigma==1000):
                R0_N_array2 = Ts_array/T_c
                ax.plot(R0_N_array2, S_SIR ,linestyle = 'dashed', linewidth = 3, color = color_p, label = r'Network')
            if(sigma==1/4):
                R0_N_array2 = Ts_array/T_c
                ax.plot(R0_N_array2, S_SEIR ,linestyle = 'dashed', linewidth = 3, color = color_p, label = r'Network')
        
    ax.hlines(1,0.5,6.5, linestyle = 'dashed', color = 'silver', alpha = .4, linewidth = 1)
    ax.vlines(1,0,1, linestyle = 'dashed', color = 'silver', alpha = .4, linewidth = 1)
    my_plot_layout(ax=ax, xlabel = r'$R_0$', ylabel=r'Probability of epidemic')
    if(sigma==1000):
        ax.set_xlim(0.7,6.05)
    if(sigma==1/4):
        ax.set_xlim(0.7,4.05)
    ax.set_ylim(-0.05, 1.05)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(np.concatenate((lines_symbols,handles)), np.concatenate((labels_symbols,labels)) , fontsize = 20, loc = 4, framealpha=.95)
    fig.savefig('../Figures/Stochastic/Networks/barabasi-albert/Prob_epi_'+labels_model[l]+'.pdf')

