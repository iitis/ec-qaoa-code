from logging import warn
import matplotlib.pyplot as plt
import numpy as np
from qaoa_xy_mixer import qaoa_xy
import matplotlib.ticker as ticker
from noise_models import create_noise_model_random_x

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif"})

##########################################################################################
##################### ###################### QAOA optimization ###########################
##########################################################################################
cities_no = [3, 4]
layers_no = 8
tsp_instances = 40
exp_no = 1
gamma = 0.002
noise_model_str = "randx"+str(gamma)
layer_init = 8
repetition_qaoa_per_layer = 3
scheme = "01"    
model = "xy"
noise_model = create_noise_model_random_x(gamma)

y_min = {3: 0.01, 4: 0.025}
y_max = {3: .06, 4: 0.08}
def opt_plot(ax, data_en_nopost, data_en_post, city_no):
    data_en_nopost = np.asarray(data_en_nopost)
    data_en_post = np.asarray(data_en_post)
    for ind, en in enumerate(data_en_post):
        if en == 0.:
            print(ind, en)

    permsort = np.asarray(np.argsort(data_en_nopost))
    warn("HACK below")
    data_en_nopost = data_en_nopost[permsort][1:]
    data_en_post = data_en_post[permsort][1:]

    assert np.max(data_en_post) < y_max[city_no], "should be {} {}".format(np.max(data_en_post), y_max[city_no])
    assert np.max(data_en_nopost) < y_max[city_no] , "should be {} {}".format(np.max(data_en_nopost), y_max[city_no])
    assert y_min[city_no] < np.max(data_en_nopost), "should be {} {}".format(y_min[city_no], np.max(data_en_nopost))
    assert y_min[city_no] < np.max(data_en_post), "should be {} {}".format(y_min[city_no], np.max(data_en_post))
    print(data_en_post/data_en_nopost)
    ax.scatter(data_en_nopost, data_en_post, c="k", marker="x", s=6, lw=0.5)
    # ax.hlines(np.min(data_en_nopost), y_min[city_no], y_max[city_no], linestyle="-", color="r", linewidth=.2)

f_nopost = dict()
f_post = dict()
normalization = "full"
f_nopost = lambda x, diag_vals, hamildict: np.real(qaoa_xy(x, noise_model, diag_vals, hamiltonian=hamildict,normalization=normalization)[0])
f_post = lambda x, diag_vals, hamildict: np.real(qaoa_xy(x, noise_model, diag_vals, midmeasure='yes', hamiltonian=hamildict, scheme=scheme, normalization=normalization)[0])

fig, axes = plt.subplots(2, 3, figsize=(6,4), sharex="row", sharey="row")
data_x = range(1, tsp_instances+1)

for city_ind, city_no in enumerate(cities_no):
    print(city_no)
    ## energy for optimized angles
    data_en_nopost = []
    data_en_post = []
    for m in range(tsp_instances):
        filename = "data_noise_qaoa_xy/optimization/tsp{}-model{}-noise{}-exp{}-scheme{}".format(city_no, model, noise_model_str, m, scheme)
        hamildict = np.load('./tsp_data_xy/tsp{}/tsp_dict_{}_{}.npz'.format(city_no, model, m+1))
        for j in [0,2]:
            for l in range(len(hamildict[:,j])):
                if hamildict[:,j][l] > 0:
                    hamildict[:,j][l] -= 1
        diag_vals = np.load('./tsp_data_xy/tsp{}/qubo_{}.npz'.format(city_no, m+1))

        angle = np.load(filename + "_nopost.npy")
        data_en_nopost.append(f_nopost(angle, diag_vals, hamildict))
        data_en_post.append(f_post(angle, diag_vals, hamildict))
    opt_plot(axes[city_ind,0], data_en_nopost, data_en_post, city_no)    

    # plt.title("Improvement of optimal solution")

    data_en_nopost = []
    data_en_post = []
    for m in range(tsp_instances):
        filename = "data_noise_qaoa_xy/optimization/tsp{}-model{}-noise{}-exp{}-scheme{}".format(city_no, model, noise_model_str, m, scheme)
        hamildict = np.load('./tsp_data_xy/tsp{}/tsp_dict_{}_{}.npz'.format(city_no, model, m+1))
        for j in [0,2]:
            for l in range(len(hamildict[:,j])):
                if hamildict[:,j][l] > 0:
                    hamildict[:,j][l] -= 1
        diag_vals = np.load('./tsp_data_xy/tsp{}/qubo_{}.npz'.format(city_no, m+1))

        angle = np.load(filename + "_nopost.npy")
        data_en_nopost.append(f_nopost(angle, diag_vals, hamildict))
        angle = np.load(filename + "_post.npy")
        data_en_post.append(f_post(angle, diag_vals, hamildict))
    opt_plot(axes[city_ind,1], data_en_nopost, data_en_post, city_no)

    data_en_nopost = []
    data_en_post = []
    for m in range(tsp_instances):
        filename = "data_noise_qaoa_xy/optimization/tsp{}-model{}-noise{}-exp{}-scheme{}".format(city_no, model, noise_model_str, m, scheme)
        hamildict = np.load('./tsp_data_xy/tsp{}/tsp_dict_{}_{}.npz'.format(city_no, model, m+1))
        for j in [0,2]:
            for l in range(len(hamildict[:,j])):
                if hamildict[:,j][l] > 0:
                    hamildict[:,j][l] -= 1
        diag_vals = np.load('./tsp_data_xy/tsp{}/qubo_{}.npz'.format(city_no, m+1))

        angle = np.load(filename + "_nopost.npy")
        data_en_nopost.append(f_nopost(angle, diag_vals, hamildict))
        angle = np.load(filename + "_nopost_post.npy")
        data_en_post.append(f_post(angle, diag_vals, hamildict))
    opt_plot(axes[city_ind,2], data_en_nopost, data_en_post, city_no)


    # axes[0].set_ylim(10e-2, 2)
axes[0,0].set_ylabel("rescaled energy mid-post\n3 cities")
axes[1,0].set_ylabel("rescaled energy mid-post\n4 cities")

axes[0,0].set_title("error mitigation")
axes[0,1].set_title("optimization")
axes[0,2].set_title("two-step opt.")

axes[1,0].set_xlabel("rescaled energy\nno-mid-post")
axes[1,1].set_xlabel("rescaled energy\nno-mid-post")
axes[1,2].set_xlabel("rescaled energy\nno-mid-post")

for ax_ind, ax in enumerate(axes.flatten()):
    ax.set_aspect('equal')

    if ax_ind < 3:
        city_no = 3
    else:
        city_no = 4
    ax.plot([y_min[city_no], y_max[city_no]], [y_min[city_no], y_max[city_no]], 'k-', alpha=0.75, zorder=0, linewidth=.2)    
    
    ax.plot([y_min[city_no], y_max[city_no]], [0.85*y_min[city_no], 0.85*y_max[city_no]], 'r-', alpha=0.75, zorder=0, linewidth=.2)    
    
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '{:.3f}'.format(y)))
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '{:.3f}'.format(y)))

    ax.set_axisbelow(True)
    for axis in [ax.xaxis, ax.yaxis]:
        axis.grid(True, which="both", linewidth=.1)
    ax.set_ylim(y_min[city_no], y_max[city_no])
    ax.set_xlim(y_min[city_no], y_max[city_no])



plt.subplots_adjust()
plt.savefig("plots/qaoa_optimization.pdf", bbox_inches="tight")
