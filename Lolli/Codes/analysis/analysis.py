#---------Import---------
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
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

betas=np.array([0.54, 0.72, 0.90])
p_in_s=np.array([0.001])
p_in_s=np.array([0.001, 0.01, 0.1])
taus=np.array([4])

colors = [plt.cm.Blues(np.linspace(0,1,int(len(betas)+1))), plt.cm.Greens(np.linspace(0,1,int(len(betas)+1))), plt.cm.Reds(np.linspace(0,1,int(len(betas)+1)))]
markers = ['^', '*', 'o']
alphas = [1, 1]
labels = [['lollipop', 'antigen', 'control'], [None, None, None]]
z = 1
z2 = 1
data_dict = {}
#fig, ax = plt.subplots(figsize = (12,10), gridspec_kw={'top':.95, 'left':.18	})
fig2, ax2 = plt.subplots(figsize = (12,10), gridspec_kw={'top':.95, 'left':.18	})
fig3, ax3 = plt.subplots(figsize = (12,10), gridspec_kw={'top':.95, 'left':.18	})
for k, p_in in enumerate(p_in_s):
	
	for i, beta in enumerate(betas):
		fig, ax = plt.subplots(figsize = (12,10), gridspec_kw={'top':.95, 'left':.18})
		for j, tau in enumerate(taus):
			means_inf = np.array([])
			means_det = np.array([])
			means_trans = np.array([])
			means_det_trans = np.array([])
			for d, n_testing_day in enumerate(n_testing_days):
				data = np.loadtxt(Text_files_path+'statistics_days-%d-'%(n_testing_day)+testing_days2[d]+'_beta-%.6f_pin-%.6f_tau-%d_%d.txt'%(beta, p_in, tau, times))
				mean_inf = ((np.mean(data[:,0])))/(1)
				mean_det = ((np.mean(data[:,1])))/(1)
				mean_trans = ((np.mean(data[:,2])))/(1)
				mean_det_trans = ((np.mean(data[:,3])))/(1)
				means_inf = np.append(means_inf, mean_inf)
				means_det = np.append(means_det, mean_det)
				means_trans = np.append(means_trans, mean_trans)
				means_det_trans = np.append(means_det_trans, mean_det_trans)
				if(n_testing_day==0):
					z = np.mean(data[:,0])
					z2 = np.mean(data[:,2])
				data_dict[(p_in, beta*5, n_testing_day)] = {'det': mean_det*100, 'det_trans': mean_det_trans*100, 'prev_inf':(z-mean_inf)/z, 'prev_trans':(z2-mean_trans)/z2}
			#--------- Infections ---------
			#--------- cumulative ---------
			#ax.plot(n_testing_days, means_inf/z,  marker = markers[k], linestyle = '--', color = colors[i][1], ms = 20, alpha = alphas[1], label = labels[0][j])
			#ax.plot(n_testing_days, means_trans/z2,  marker = markers[k], linestyle = '-', color = colors[i][2], ms = 20, alpha = alphas[0])
			
			#my_plot_layout(ax=ax, yscale = 'linear', xscale = 'linear', ticks_labelsize = 28,
	        #           xlabel = 'Tests per week', ylabel = 'Cumulative infections', title = '', x_fontsize=28, y_fontsize = 28,
	        #           t_fontsize = 24)
			#ax.set_xticks(n_testing_days)
			#ax.yaxis.set_major_formatter(PercentFormatter(1))
			#ax.legend(fontsize = 24, loc = 0)

			#--------- Infections ---------
			#--------- prevented ---------
			ax2.plot(n_testing_days, (z-means_inf)/z,  marker = markers[k], linestyle = '--', color = colors[i][1], ms = 20, alpha = alphas[1], label = labels[0][j])
			ax2.plot(n_testing_days, (z2-means_trans)/z2,  marker = markers[k], linestyle = '-', color = colors[i][2], ms = 20, alpha = alphas[0])
			my_plot_layout(ax=ax2, yscale = 'linear', xscale = 'linear', ticks_labelsize = 28,
	                   xlabel = 'Tests per week', ylabel = 'Prevented infections', title = '', x_fontsize=28, y_fontsize = 28,
	                   t_fontsize = 24)
			ax2.set_xticks(n_testing_days)
			ax2.yaxis.set_major_formatter(PercentFormatter(1))
			ax2.legend(fontsize = 24, loc = 0)
			
			#--------- Detections ---------
			ax3.plot(n_testing_days, (means_det)*100,  marker = markers[k], linestyle = '--', color = colors[i][2], ms = 20, alpha = alphas[0], label = labels[0][j])
			my_plot_layout(ax=ax3, yscale = 'log', xscale = 'linear', ticks_labelsize = 28,
	                   xlabel = 'Tests per week', ylabel = 'Detections', title = '', x_fontsize=28, y_fontsize = 28,
	                   t_fontsize = 24)
			ax3.set_xticks(n_testing_days)
			ax3.legend(fontsize = 24, loc = 0)

			ax.bar(n_testing_days, (means_det-means_det_trans)*100, color = colors[i][2], label = 'Imported')
			ax.bar(n_testing_days, (means_det_trans)*100, color = colors[i][1], bottom=(means_det-means_det_trans)*100, label = 'Trans.')
			my_plot_layout(ax=ax, yscale = 'linear', xscale = 'linear', ticks_labelsize = 28,
	                   xlabel = 'Tests per week', ylabel = 'Detections', title = r'$p_{in} = %.1f\%%$ ; $R_0 = %.1f$'%(p_in*100, beta*5), x_fontsize=28, y_fontsize = 28,
	                   t_fontsize = 24)
			ax.legend(loc = 2, fontsize=24)
			fig.savefig('../../Figures/Kitas_Schools/Analysis/detections_pin-%.3f_beta-%.3f.png'%(p_in, beta))

		#fig.savefig('../../Figures/Kitas_Schools/newplots/statistics_infections_percent_pin-%.3f_beta-%.3f.png'%(p_in, beta))
		#fig2.savefig('../../Figures/Kitas_Schools/newplots/statistics_infections_prevented_percent_pin-%.3f_beta-%.3f.png'%(p_in, beta))
		#fig3.savefig('../../Figures/Kitas_Schools/newplots/statistics_detections_pin-%.3f_beta-%.3f.png'%(p_in, beta))

data_frame = pd.DataFrame.from_dict(data_dict)
data_frame.to_csv(Text_files_path+'data_simulations.csv')
#data_frame.to_excel(Text_files_path+'data_simulations.xlsx')

#--------- Legend 1 ---------
custom_lines = [Line2D([0], [0],linestyle = '', marker = 's', ms = 20, color=colors[i][2], lw=4) for i in range(int(len(betas)))]
labels = [r'$%.1f$'%(beta*5) for beta in betas]

#legend1_1 = ax.legend(custom_lines, labels, fontsize = 26, loc = 1, title = r'$p_{in}$', title_fontsize = 24, framealpha = .6)
#legend1_2 = ax2.legend(custom_lines, labels, fontsize = 26, loc = 2, title = r'$p_{in}$', title_fontsize = 24, framealpha = .6)
legend1_3 = ax3.legend(custom_lines, labels, fontsize = 26, loc = 2, title = r'$R_0$', title_fontsize = 24, framealpha = .6)

#ax.add_artist(legend1_1)
#ax2.add_artist(legend1_2)
ax3.add_artist(legend1_3)

#----------------------------

#--------- Legend 2 ---------
custom_lines = [Line2D([0], [0],linestyle = '--', marker = 's', ms = 10, color='k', lw=4),  Line2D([0], [0],linestyle = '-', marker = 's', ms = 10, color='k', lw=4)]
labels = ['Total', 'Transmisions']

legend2_1 = ax.legend(custom_lines, labels, fontsize = 26, loc = 1, title = 'Type', title_fontsize = 24, framealpha = .6)
legend2_2 = ax2.legend(custom_lines, labels, fontsize = 26, loc = 4, title = 'Type', title_fontsize = 24, framealpha = .6)
#legend2_3 = ax3.legend(custom_lines, labels, fontsize = 24, loc = 2, title = 'Type', title_fontsize = 24)

ax.add_artist(legend2_1)
ax2.add_artist(legend2_2)
#ax3.add_artist(legend2_3)
#----------------------------
#--------- Legend 3 ---------
custom_lines = [Line2D([0], [0],linestyle = '', marker = i, ms = 20, color='k', lw=4) for i in np.flip(markers)]
labels = [r'$%.1f\%%$'%(p_in*100) for p_in in np.flip(p_in_s)]

#ax.legend(custom_lines, labels, fontsize = 24, loc = 3, title = r'$R_0$', title_fontsize = 24, framealpha = .6)
#ax2.legend(custom_lines, labels, fontsize = 24, loc = 4, title = r'$R_0$', title_fontsize = 24, framealpha = .6)
ax3.legend(custom_lines, labels, fontsize = 26, loc = 3, title = r'$p_{in}$', title_fontsize = 24, framealpha = .6)

#fig.savefig('../../Figures/Kitas_Schools/Analysis/statistics_infections_percent.png')
fig2.savefig('../../Figures/Kitas_Schools/Analysis/statistics_infections_prevented_percent.png')
fig3.savefig('../../Figures/Kitas_Schools/Analysis/statistics_detections.png')

