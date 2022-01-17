#%%
from math import log
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
import sys, os
import numpy as np

mpl.rc('font', family='serif')
mpl.rc('text', usetex=True)

if len(sys.argv) < 2:
    print("You should put some inputs here. If you don't know, type '-h','--help' for instructions.")
    exit(0)

if sys.argv[1] in ["--help", "-h"]:
    print()
    print(" -> Please specify:")
    print()
    print("   - number of: layers, instances, experiments;")
    print()
    print("- deviation calculation ('rand or optim')")
    print()
    print("Example: python energy_diff_plot.py 20 10 10 optim")
    exit(0)


model = "xy"

# noise_model_str = sys.argv[2]
# assert noise_model_str in ['amp_damp','depol','rand_x'], "Noise model should be 'amp_damp', 'depol', 'rand_x'."
# if noise_model_str == 'amp_damp':
#     noise_md_st = "aplitude damping"
# if noise_model_str == 'depol':
#     noise_md_st = "depolarazing channel"
# if noise_model_str == 'rand_x':
#     noise_md_st = "random unitary x"

# Fix naming
# city_no = int(sys.argv[2])
layers_no = int(sys.argv[1]) # 10 int(len(angles_list)/2)
instances_no = int(sys.argv[2]) # 5 or 10exp_no = 1
exp_no = int(sys.argv[3])

angles_typ = sys.argv[4]
assert angles_typ in ['rand','optim'], "angles_typ should be 'rand' or 'optim'"
if angles_typ == 'rand':
    angls_typ = 'random'
if angles_typ == 'optim':
    angls_typ = 'optimized'

mean_calc = 'perc'


print(" - The model is",model)
print(" - cities,",layers_no,"layers,",instances_no ,"intances," "and",exp_no,"experiments")
print(" - Angles type:", angls_typ)
print()
#%%
## loading data
data_prob = {}
data_prob = {}
noise_model_str = ['amp_damp', 'depol', 'rand_x']
noise_lebel = ['amp. damping', 'depolarizing', 'bit flip']
# diag_vals = 
my_keys = [None, 'yes']
#noises = [0.01, 0.005, 0.002, 0.001, 0.0005]
noises = [0.01, 0.005, 0.002, 0.0005]
fig, ax = plt.subplots(2, 3, sharex=True, sharey="row")
list_types = ['b-', 'g-.', 'm:', 'r--', 'c-']
list_color = [el[0] for el in list_types]
percentile = 10
assert 0 <= percentile <= 50
city_no = [3,4]
scheme = "0001"
for cit, city in enumerate(city_no):
    for m, noise_model in enumerate(noise_model_str):
        for p in noises:
            for i in range(1,exp_no+1):
                data_prob[p,i] = np.load('data_noise_qaoa_xy/tsp{}_plotting_for_paper/{}-model-{}-noise-{}-exp_{}-scheme{}-angles-keys_layers{}_instances-prob.npy'.format(city, model,
                                                                                                                                            noise_model, p, i, scheme, layers_no))
        # #%%
        ## calculating mean (over angles) of relative energy/prob
        data_prob_mean = {}
        for p in noises:
            for i in range(1,exp_no+1):
                for k_ind in range(len(my_keys)):
                    data_prob_mean[p,i,k_ind] = np.zeros(layers_no)
                    for a in range(instances_no):
                        nominator = data_prob[p,i][k_ind+1,:,a]
                        denominator = data_prob[p,i][-1,:,a]
                        data_prob_mean[p,i,k_ind] += denominator
                    data_prob_mean[p,i,k_ind] /= instances_no

        #
        # calculating statistics over TSP instances
        data_prob_mean_plot = {}
        data_prob_std_plot = {}
        data_prob_min_plot = {}
        data_prob_max_plot = {}

        for p in noises:
            for k_ind in range(len(my_keys)):
                if mean_calc == 'std':
                    data_prob_mean_plot[p,k_ind] = sum([data_prob_mean[p,i,k_ind] for i in range(1,exp_no+1)]) / exp_no
                    data_prob_std_plot[p,k_ind] = np.sqrt(sum([(data_prob_mean[p,i,k_ind] - data_prob_mean_plot[p,k_ind])**2 for i in range(1,exp_no+1)])) / np.sqrt(exp_no)
                    data_prob_min_plot[p,k_ind] = data_prob_mean_plot[p,k_ind] - data_prob_std_plot[p,k_ind]
                    data_prob_max_plot[p,k_ind] = data_prob_mean_plot[p,k_ind] + data_prob_std_plot[p,k_ind]
                elif mean_calc == 'ex':
                    data_prob_mean_plot[p,k_ind] = sum([data_prob_mean[p,i,k_ind] for i in range(1,exp_no+1)]) / exp_no
                    data_prob_max_plot[p,k_ind] = [np.max([data_prob_mean[p,i,k_ind][t] for i in range(1,exp_no+1)]) for t in range(layers_no)]
                    data_prob_min_plot[p,k_ind] = [np.min([data_prob_mean[p,i,k_ind][t] for i in range(1,exp_no+1)]) for t in range(layers_no)]
                elif mean_calc == 'perc':
                    data_prob_mean_plot[p,k_ind] = sum([data_prob_mean[p,i,k_ind] for i in range(1,exp_no+1)]) / exp_no
                    data_prob_max_plot[p,k_ind] = [np.percentile([data_prob_mean[p,i,k_ind][t] for i in range(1,exp_no+1)], 100-percentile) for t in range(layers_no)]
                    data_prob_min_plot[p,k_ind] = [np.percentile([data_prob_mean[p,i,k_ind][t] for i in range(1,exp_no+1)], percentile) for t in range(layers_no)]
                else:
                    raise ValueError(f"Incorrect value for mean_calc {mean_calc}")
                
    #
    #with filling plotting!!
    # y_prob = {}
    # y_prob_up = {}
    # y_prob_down = {}
        x = range(1,layers_no+1)
        y = [0]*layers_no
        skip = 3
        for j in range(len(noises)):
            y_prob = [data_prob_mean_plot[noises[j], 0]][0]
            y_prob_up = [data_prob_max_plot[noises[j],0]][0]
            y_prob_down = [data_prob_min_plot[noises[j],0]][0]
            ax[cit,m].plot(x[::skip],y_prob[::skip], list_types[j])
            ax[cit,m].fill_between(x[::skip], y_prob_up[::skip], y_prob_down[::skip], alpha=0.4, color = list_color[j])
            # ax[cit,m].set_yscale('log')

        ax[0,0].text(0.1, 0.1, '(a)', ha='center', va='center', transform=ax[0,0].transAxes)
        ax[0,1].text(0.1, 0.1, '(b)', ha='center', va='center', transform=ax[0,1].transAxes)
        ax[0,2].text(0.1, 0.1, '(c)', ha='center', va='center', transform=ax[0,2].transAxes)
        ax[1,0].text(0.1, 0.1, '(d)', ha='center', va='center', transform=ax[1,0].transAxes)
        ax[1,1].text(0.1, 0.1, '(e)', ha='center', va='center', transform=ax[1,1].transAxes)
        ax[1,2].text(0.1, 0.1, '(f)', ha='center', va='center', transform=ax[1,2].transAxes)
        
        # for a in ax[0,:]:
        #     a.set_ylim(0, 1)
        #     a.set_yticks([])
        # ax[0,0].set_yticks([0,0.01,0.02,0.03])

        # for a in ax[1,:]:
        #     a.set_ylim(0, 1)
        #     a.set_yticks([])
        # ax[1,0].set_yticks([-0.001,0,0.001,0.002,0.003])



        # ax[1,0].set_ylabel('$\overline{P}$', fontsize=15)

for nse, noise in enumerate(noise_lebel):
    ax[0,nse].set_title('{}'.format(noise), fontsize = 10)

ax[0,0].set_yscale("log")
ax[1,0].set_yscale("log")
ax[0,0].set_ylim(1e-5,3)
ax[1,0].set_ylim(5e-11,3)
ax[0,0].set_yticks([1,1e-1,1e-2,1e-3,1e-4])

for cit in range(len(city_no)):
    for m in range(len(noise_model_str)):
        ax[cit,m].plot([1,40], [1,1], 'k--', linewidth=.7)

# fig.show()
# plt.yticks([0,0.5*1e-1,1e-1, 1e-2])


fig.legend(title = 'Noise $(\gamma) = $' , loc='upper center', labels = noises,  borderaxespad=2.5, ncol=5)
fig.text(0.5, 0.02, 'Number of levels', ha='center', fontsize = 12)
fig.text(0.02, 0.5, 'Probability of accepting a circuits run', va='center', rotation='vertical', fontsize = 12)
# fig.tight_layout()
fig.subplots_adjust(top=0.76)
fig.savefig("plots/{}_city_{}-model-{}-noise--angles_plot_prob.pdf".format(city_no, model, angles_typ))

