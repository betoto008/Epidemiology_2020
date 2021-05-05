#---------Import---------
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib.ticker import PercentFormatter
from matplotlib.lines import Line2D
#---------Functions---------
def my_linear_func(x, a, b):
    return a + b*x
def my_quadratic_func(x, a, b, c):
    return np.log(a)+np.log(np.sqrt(-b)) + b*(x-c)**2
def my_plot_layout(ax, yscale = 'linear', xscale = 'linear', ticks_labelsize = 24,
                   xlabel = '', ylabel = '', title = '', x_fontsize=24, y_fontsize = 24,
                   t_fontsize = 24):
    ax.tick_params(labelsize = ticks_labelsize)
    ax.set_yscale(yscale)
    ax.set_xscale(xscale)
    ax.set_xlabel(xlabel, fontsize = x_fontsize)
    ax.set_ylabel(ylabel, fontsize = y_fontsize)
    ax.set_title(title, fontsize = y_fontsize)
#------------------------------------

Text_files_path = '../../../../../Dropbox/Research/Epidemiology_2020/Text_files/Kitas_Schools/Statistics/'

days = ["mon", "tue", "wed", "thu", "fri", "sat", "sun"]
T = 7*4
times = 1

what = ['cum_inf', 'cum_det', 'prev_inf']

n_testing_days=[0, 1, 2, 3]
testing_days=[[0], [2], [1,3], [0, 2, 4]]
testing_days2=["0", "2", "13", "024"]

betas=[0.36, 0.54, 0.72]
p_in_s=[0.0011, 0.011, 0.11]
taus=[4]

colors = [plt.cm.Purples(np.linspace(0,1,int(len(betas)+1))), plt.cm.Blues(np.linspace(0,1,int(len(betas)+1))), plt.cm.Greens(np.linspace(0,1,int(len(betas)+1))), plt.cm.Reds(np.linspace(0,1,int(len(betas)+1)))]
markers = ['^', '*', 'o', 's', '.']
alphas = [1, 0.5]
labels = [['lollipop', 'antigen', 'control'], [None, None, None]]
z = 1
z2 = 1

fig, ax = plt.subplots(figsize = (12,10), gridspec_kw={'top':.95, 'left':.18	})
fig2, ax2 = plt.subplots(figsize = (12,10), gridspec_kw={'top':.95, 'left':.18	})
fig3, ax3 = plt.subplots(figsize = (12,10), gridspec_kw={'top':.95, 'left':.18	})
for k, p_in in enumerate(p_in_s):
	for i, beta in enumerate(betas):
		for j, tau in enumerate(taus):
			means_inf = np.array([])
			means_det = np.array([])
			means_trans = np.array([])
			for d, n_testing_day in enumerate(n_testing_days):
				data = np.loadtxt(Text_files_path+'statistics_days-%d-'%(n_testing_day)+testing_days2[d]+'_beta-%.6f_pin-%.6f_tau-%d_%d.txt'%(beta, p_in, tau, times))
				means_inf = np.append(means_inf, ((np.mean(data[:,0])))/(1))
				means_det = np.append(means_det, ((np.mean(data[:,1])))/(1))
				means_trans = np.append(means_trans, ((np.mean(data[:,2])))/(1))
				if(n_testing_day==0):
					z = np.mean(data[:,0])
					z2 = np.mean(data[:,2])
			#--------- Infections ---------
			#--------- cumulative ---------
			plot = ax.plot(n_testing_days, means_inf/z,  marker = markers[0], linestyle = '--', color = colors[k][i+1], ms = 20, alpha = alphas[0], label = labels[0][j])
			ax.plot(n_testing_days, means_trans/z2,  marker = markers[1], linestyle = '--', color = colors[k][i+1], ms = 20, alpha = alphas[0])
			
			my_plot_layout(ax=ax, yscale = 'linear', xscale = 'linear', ticks_labelsize = 28,
	                   xlabel = 'Tests per week', ylabel = 'Cumulative infections', title = '', x_fontsize=28, y_fontsize = 28,
	                   t_fontsize = 24)
			ax.set_xticks(n_testing_days)
			ax.yaxis.set_major_formatter(PercentFormatter(1))
			ax.legend(fontsize = 24, loc = 0)

			#--------- Infections ---------
			#--------- prevented ---------
			ax2.plot(n_testing_days, (z-means_inf)/z,  marker = markers[0], linestyle = '--', color = colors[k][i+1], ms = 20, alpha = alphas[0], label = labels[0][j])
			ax2.plot(n_testing_days, (z2-means_trans)/z2,  marker = markers[1], linestyle = '--', color = colors[k][i+1], ms = 20, alpha = alphas[0])
			my_plot_layout(ax=ax2, yscale = 'linear', xscale = 'linear', ticks_labelsize = 28,
	                   xlabel = 'Tests per week', ylabel = 'Prevented infections', title = '', x_fontsize=28, y_fontsize = 28,
	                   t_fontsize = 24)
			ax2.set_xticks(n_testing_days)
			ax2.yaxis.set_major_formatter(PercentFormatter(1))
			ax2.legend(fontsize = 24, loc = 0)
			
			#--------- Detections ---------
			ax3.plot(n_testing_days, (means_det)*100,  marker = markers[2], linestyle = '--', color = colors[k][i+1], ms = 20, alpha = alphas[0], label = labels[0][j])
			my_plot_layout(ax=ax3, yscale = 'log', xscale = 'linear', ticks_labelsize = 28,
	                   xlabel = 'Tests per week', ylabel = 'Detections', title = '', x_fontsize=28, y_fontsize = 28,
	                   t_fontsize = 24)
			ax3.set_xticks(n_testing_days)
			ax3.legend(fontsize = 24, loc = 0)


		#fig.savefig('../../Figures/Kitas_Schools/newplots/statistics_infections_percent_pin-%.3f_beta-%.3f.png'%(p_in, beta))
		#fig2.savefig('../../Figures/Kitas_Schools/newplots/statistics_infections_prevented_percent_pin-%.3f_beta-%.3f.png'%(p_in, beta))
		#fig3.savefig('../../Figures/Kitas_Schools/newplots/statistics_detections_pin-%.3f_beta-%.3f.png'%(p_in, beta))

custom_lines = [Line2D([0], [0],linestyle = '', marker = 'o', ms = 20, color=colors[i][2], lw=4) for i in range(int(len(p_in_s)))]

labels = [r'$%.1f\%%$'%(p_in*100) for p_in in p_in_s]

legend1_1 = ax.legend(custom_lines, labels, fontsize = 24, loc = 1, title = r'$p_{in}$', title_fontsize = 24)
legend1_2 = ax2.legend(custom_lines, labels, fontsize = 24, loc = 1, title = r'$p_{in}$', title_fontsize = 24)
legend1_3 = ax3.legend(custom_lines, labels, fontsize = 24, loc = 1, title = r'$p_{in}$', title_fontsize = 24)

ax.add_artist(legend1_1)
ax2.add_artist(legend1_2)
ax3.add_artist(legend1_3)

custom_lines = [Line2D([0], [0],linestyle = '', marker = 's', ms = 20, color=i, lw=4) for i in plt.cm.Greys(np.linspace(0,1,int(len(betas)+1))[1:])]

labels = [r'$%.1f$'%(beta*5) for beta in betas]

ax.legend(custom_lines, labels, fontsize = 24, loc = 3, title = r'$R_0$', title_fontsize = 24)
ax2.legend(custom_lines, labels, fontsize = 24, loc = 3, title = r'$R_0$', title_fontsize = 24)
ax3.legend(custom_lines, labels, fontsize = 24, loc = 3, title = r'$R_0$', title_fontsize = 24)

fig.savefig('../../Figures/Kitas_Schools/Analysis/statistics_infections_percent.png')
fig2.savefig('../../Figures/Kitas_Schools/Analysis/statistics_infections_prevented_percent.png')
fig3.savefig('../../Figures/Kitas_Schools/Analysis/statistics_detections.png')

