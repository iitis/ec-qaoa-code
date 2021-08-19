import multiprocessing
import os
import sys
from time import time

from scipy.optimize import minimize

from noise_models import *
from qaoa_xy_mixer import qaoa_xy

city_no = int(sys.argv[1])
layers_no = 8
tsp_instances = int(sys.argv[2])

gamma = 0.002
noise_model_str = "randx"+str(gamma)
layer_init = layers_no

model = "xy"
scheme = "01"

pool_size = 20
noise_model = create_noise_model_random_x(gamma)

#noise strength

def f(m):
    filename = "data_noise_qaoa_xy/optimization/tsp{}-model{}-noise{}-exp{}-scheme{}".format(city_no, model, noise_model_str, m, scheme)
    if os.path.exists(filename + "_nopost.npy"):
        print(filename)
        print("({}) > File exists - skipped".format(filename))
        return True

    hamildict = np.load('./tsp_data_xy/tsp{}/tsp_dict_{}_{}.npz'.format(city_no, model, m+1))
    for j in [0,2]:
        for l in range(len(hamildict[:,j])):
            if hamildict[:,j][l] > 0:
                hamildict[:,j][l] -= 1
    diag_vals = np.load('./tsp_data_xy/tsp{}/qubo_{}.npz'.format(city_no, m+1))
    f_nopost = lambda x: np.real(qaoa_xy(x, noise_model, diag_vals, hamiltonian=hamildict, normalization="shift")[0])
    f_post = lambda x: np.real(qaoa_xy(x, noise_model, diag_vals, midmeasure='yes', hamiltonian=hamildict, scheme=scheme, normalization="shift")[0])

    np.random.seed()
    
    angles = np.random.uniform(high=2*np.pi, size=2*layer_init)
    minimizer = lambda f, x: minimize(f, x, method="L-BFGS-B")
    print("Experiment {} started".format(m))
    
    def trajectory_optimizer(angle, f):
        result = minimizer(f, angle)
        for layers in range(layer_init+1, layers_no):
            angle = result.x
            qubo_theta = angles[0:layers]
            mixer_phi = angles[layers:]
            angle = np.concatenate([qubo_theta, [np.random.uniform(0, 2*np.pi)], mixer_phi, [np.random.uniform(0, 2*np.pi)]])

            result = minimizer(f, angle)
        return result
    
    print("Optimization started ({}, {})".format(m, filename))

    t1 = time()
    res_nopost = trajectory_optimizer(angles, f_nopost)
    print("No-postselection done ({}, {}, {})".format(m, res_nopost.fun, time()-t1))
    print(res_nopost.x / np.pi)

    t1 = time()
    res_post = trajectory_optimizer(angles, f_post)
    print("Postselection done ({}, {}, {})".format(m, res_post.fun - res_nopost.fun, time()-t1))
    print(res_post.x / np.pi)

    t1 = time()
    res_post_nopost = minimizer(f_nopost, res_post.x)
    print("res_post_nopost done ({}, {}, {})".format(m, res_post_nopost.fun - res_post.fun, time()-t1))

    t1 = time()
    res_nopost_post = minimizer(f_post, res_nopost.x)
    print("res_nopost_post done ({}, {}, {})".format(m, res_nopost_post.fun - res_nopost.fun, time()-t1))
    
    np.save(filename + "_nopost", res_nopost.x)
    np.save(filename + "_post", res_post.x)
    np.save(filename + "_post_nopost", res_post_nopost.x)
    np.save(filename + "_nopost_post", res_nopost_post.x)
    print("############################ {} done #####################".format(m))
    return True

pool = multiprocessing.Pool(pool_size)
pool.map(f, range(tsp_instances))
