import matplotlib.pyplot as plt
from encoding_qubo_decoding_mixture import *
import sys, os
import random

city_no = 3

def create_noise_model(p):
    error = noise.amplitude_damping_error(p)
    error2 = error.tensor(error)
    noise_model = noise.NoiseModel()
    noise_model.add_all_qubit_quantum_error(error, ['u1', 'u2', 'u3'])
    noise_model.add_all_qubit_quantum_error(error2, ['cx'])
    return noise_model

## Energy Obtaining Parameters ##
times = 30 # 10 int(len(angles_list)/2)
no_angles = 50 # 20 or more
instances_no = 50 # 5 or 10

#my_keys = [(None, False), ('end', False), ('end', 'no-ps'), ('end', 'ps')]
my_keys = [(None, False), ('end', False), ('end', 'ps')]
print(sys.argv[1])

p = [0.005, 0.01, 0.02][int(sys.argv[1])]

for m in range(instances_no):
    print("({}) m = {}".format(p, m))
    filename = "data_noise/amp_damp_{}-inst_{}-keys_angles_times-".format(p, m)
    if all(os.path.exists(filename + end + ".npy") for end in ["en", "prob"]):
        print("({}) > File exists - skipped")
        continue

    hamildict = np.load('./qaoa_efficiency_analysis/tsp_dict/tsp_{}-{}.npz'.format(city_no, m+1))
    diag_vals = np.load('./qaoa_efficiency_analysis/tsp_dict/qubo_vec_{}-{}.npz'.format(city_no, m+1))
    
    data_en = np.zeros((len(my_keys)+1, no_angles, times))
    data_prob = np.zeros((len(my_keys)+1, no_angles, times))
    for a in range(no_angles):
        print("({}) > {}-th angle".format(p, a))
        angles_list = np.random.random(2*times) * 2 * np.pi

        #outputs energy
        f_pure = lambda t: hobo_qubo_mixer_measurement_pure(hamildict, angles_list[0:2*t], diag_vals)[0]
        data_en[0,a,:] = np.asmatrix([f_pure(t) for t in range(times)])
        data_en[1,a,:] = np.ones(times)
        
        f = lambda t, x, y: hobo_qubo_mixer_measurement(hamildict, create_noise_model(p), angles_list[0:2*t], diag_vals, x, y)
        for (k_ind, k) in enumerate(my_keys):
            data_raw = [f(t, *k) for t in range(times)]
            data_en[k_ind+1,a,:] = np.asmatrix([el[0] for el in data_raw])
            data_prob[k_ind+1,a,:] = np.asmatrix([el[1] for el in data_raw])
    
    np.save(filename + "en", data_en)
    np.save(filename + "prob", data_prob)
