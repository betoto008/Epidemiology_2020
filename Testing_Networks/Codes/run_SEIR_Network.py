from models import *
from functions import *
import matplotlib.animation as animation
from celluloid import Camera


N = 50
emp = 40
pat = N-emp
constructor = [(emp, N*5, (1)/(2)), (pat, pat, 1/5)]
G = networkx.random_shell_graph(constructor)

model = SEIRSNetworkModel(G = G, beta = 1/2, sigma = 1, gamma = 1/6, p = 0, initI = 1, store_Xseries = True)
Xseries, Tseries = model.run(T=300)

print(len(Xseries[:,0]))
choices = ['g', 'y', 'r', 'b']

fig, ax = plt.subplots(figsize=(12,8))

#camera = Camera(fig)

#for i in range(len(Xseries[:,0])):
#	conditions = [Xseries[i,:]==1,Xseries[i,:]==2,Xseries[i,:]==3,Xseries[i,:]==6]
#	networkx.draw_shell(G, nlist = [range(emp), range(emp, pat+emp)], node_color=np.select(conditions,choices,'k'), style ='dashed', ax=ax)
#	plt.pause(0.01)
#	camera.snap()
def animate(i):
    conditions = [Xseries[i,:]==1,Xseries[i,:]==2,Xseries[i,:]==3,Xseries[i,:]==6]
    networkx.draw_shell(G, nlist = [range(emp), range(emp, pat+emp)], node_color=np.select(conditions,choices,0), style ='dashed', width = 0.5, node_size = 50, ax=ax)

ani = animation.FuncAnimation(fig, animate, frames=len(Xseries[:,0]), repeat = False)

#ani = camera.animate()
ani.save('../Figures/animation.gif', writer='pillow', fps=4, dpi = 50)

