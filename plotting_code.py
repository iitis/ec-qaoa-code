#%%
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})
## for Palatino and other serif fonts use:

mpl.rc('text', usetex = True)

times = 30 # 10 int(len(angles_list)/2)
no_angles = 50 # 20 or more
instances_no = 50 # 5 or 10
noises = [0.005, 0.01, 0.02]

#%%
## loading data
data_en = {}
data_prob = {}
my_keys = [(None, False), ('end', False), ('end', 'ps')] # the last one is the relative one
for p in noises:
    for i in range(instances_no):
        data_en[p,i] = np.load('data_noise/amp_damp_{}-inst_{}-keys_angles_times-en.npy'.format(p, i))
        data_prob[p,i] = np.load('data_noise/amp_damp_{}-inst_{}-keys_angles_times-prob.npy'.format(p, i))
        

#%%
## calculating mean (over angles) of relative energy/prob
data_en_mean = {}
data_prob_mean = {}
for p in noises:
    for i in range(instances_no):
        for k_ind in range(len(my_keys)):
            data_en_mean[p,i,k_ind] = np.zeros(times)
            data_prob_mean[p,i,k_ind] = np.zeros(times) 
            for a in range(no_angles):
                nominator = np.abs(data_en[p,i][k_ind+1,a,:] - data_en[p,i][0,a,:])
                denominator = np.abs(data_en[p,i][-1,a,:] - data_en[p,i][0,a,:])
                data_en_mean[p,i,k_ind] += nominator/denominator

                data_prob_mean[p,i,k_ind] += data_prob[p,i][k_ind+1,a,:]
            data_en_mean[p,i,k_ind] /= no_angles
            data_prob_mean[p,i,k_ind] /= no_angles
#%%
## calculating statistics over TSP instances
data_en_mean_plot = {}
data_en_std_plot = {}
data_en_min_plot = {}
data_en_max_plot = {}

data_prob_mean_plot = {}
data_prob_std_plot = {}
data_prob_min_plot = {}
data_prob_max_plot = {}
for p in noises:
    for k_ind in range(len(my_keys)):
        data_en_mean_plot[p,k_ind] = sum([data_en_mean[p,i,k_ind] for i in range(instances_no)]) / instances_no
        data_en_std_plot[p,k_ind] = np.sqrt(sum([(data_en_mean[p,i,k_ind] - data_en_mean_plot[p,k_ind])**2 for i in range(instances_no)])) / np.sqrt(instances_no)
        
        data_prob_mean_plot[p,k_ind] = sum([data_prob_mean[p,i,k_ind] for i in range(instances_no)]) / instances_no
        data_prob_std_plot[p,k_ind] = np.sqrt(sum([(data_prob_mean[p,i,k_ind] - data_prob_mean_plot[p,k_ind])**2 for i in range(instances_no)])) / np.sqrt(instances_no)

        #data_en_min_plot[p,k_ind] = data_en_mean_plot[p,k_ind] - data_en_std_plot[p,k_ind] 
        #data_en_max_plot[p,k_ind] = data_en_mean_plot[p,k_ind] + data_en_std_plot[p,k_ind]

        data_en_max_plot[p,k_ind] = [np.max([data_en_mean[p,i,k_ind][t] for i in range(instances_no)]) for t in range(times)]
        data_en_min_plot[p,k_ind] = [np.min([data_en_mean[p,i,k_ind][t] for i in range(instances_no)]) for t in range(times)]

        #data_prob_min_plot[p,k_ind] = data_prob_mean_plot[p,k_ind] - data_prob_std_plot[p,k_ind] 
        #data_prob_max_plot[p,k_ind] = data_prob_mean_plot[p,k_ind] + data_prob_std_plot[p,k_ind] 

        data_prob_max_plot[p,k_ind] = [np.max([data_prob_mean[p,i,k_ind][t] for i in range(instances_no)]) for t in range(times)]
        data_prob_min_plot[p,k_ind] = [np.min([data_prob_mean[p,i,k_ind][t] for i in range(instances_no)]) for t in range(times)]


# #%%
# #with filling plotting!!
list_kraus = ['no postselection', 'postselection in the end', 'postselection in the end and inbetween']
list_color = ['r', 'b', 'g']
list_types = ['r--', 'b-', 'g:']

fig = mpl.figure.Figure(figsize=(8,4.5))
ax = fig.subplots(2,len(noises))

x = range(1,times+1)
for j in range(len(noises)):
    for n in range(len(my_keys)):
        y_en = [data_en_mean_plot[noises[j],n]][0]
        y_en_up = [data_en_max_plot[noises[j],n]][0]
        y_en_down = [data_en_min_plot[noises[j],n]][0]
        ax[0,j].plot(x,y_en, list_types[n], label = list_kraus[n])
        ax[0,j].set_title('$\gamma= {}$'.format(noises[j]))
        ax[0,j].fill_between(x, y_en_up, y_en_down, color=list_color[n], alpha=0.2)
        ax[0,j].set_yscale('log',base=10)
        ax[0,j].set_ylim(0.8,9000)
        ax[0,j].set_xlim(-0.5, np.max(x)+1.5)
        ax[0,j].hlines([10,100,1000,10_000], -0.5, times+1.5, linestyle=":", linewidth=.3, color="gray")
        # ax[0,j].hlines(1, 1, times+1, linestyle="--", color="k")
    for n in range(len(my_keys)):
        y_prob = [data_prob_mean_plot[noises[j],n]][0]
        y_prob_up = [data_prob_max_plot[noises[j],n]][0]
        y_prob_down = [data_prob_min_plot[noises[j],n]][0]
        ax[1,j].plot(x,y_prob, list_types[n])
        ax[1,j].set_xlabel('number of levels', fontsize=11)
        ax[1,j].set_yscale('log',base=10)
        #ax[1,j].fill_between(x, y_prob_up, y_prob_down, color=list_color[n], alpha=0.2)
        ax[1,j].set_xlim(-0.5, np.max(x)+1.5)
        ax[1,j].set_ylim(1e-6, 2)


ax[0,0].set_ylabel('$\Delta E$', fontsize=11)
ax[1,0].set_ylabel('$P$', fontsize=11)
    
fig.legend(loc='upper center', borderaxespad=0.5, labels=list_kraus, ncol=3)
fig.tight_layout()
fig.subplots_adjust(top=0.87)    
fig.savefig("plots/noise_plot.pdf")