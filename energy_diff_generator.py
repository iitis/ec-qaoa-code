from qaoa_xy_mixer import *
import sys, os
from noise_models import *

if len(sys.argv) < 2:
    print("You should put some inputs here. If you don't know, type '-h','--help' for instructions.")
    exit(0)

if sys.argv[1] in ["--help", "-h"]:
    print()
    print(" -> Please specify:")
    print()
    print("   - noise model ('amp_damp', 'depol' or 'rand_x');")
    print("   - number of: cities, layers, intances, experiments;")
    print("   - and type of angles ('rand' or 'optim').")
    print("   - and noise strength number (0-4)")
    print()
    print("Exemple: python energy_generator_xy.py xy amp_damp 3 20 10 10 optim")
    exit(0)

model = "xy"
noise_model_str = sys.argv[1]
assert noise_model_str in ['amp_damp','depol','rand_x'], "Noise model should be 'amp_damp', 'depol', 'rand_x'."
if noise_model_str == 'amp_damp':
    noise_md_st = "amplitude damping"
    nse_md = create_noise_model_amp_damping
if noise_model_str == 'depol':
    noise_md_st = "depolarizing channel"
    nse_md = create_noise_model_depolarazing
if noise_model_str == 'rand_x':
    noise_md_st = "random unitary x"
    nse_md = create_noise_model_random_x

city_no = int(sys.argv[2])

# ## Energy Obtaining Parameters ##
layers_no = int(sys.argv[3]) # 10 int(len(angles_list)/2)
instances_no = int(sys.argv[4]) # 5 or 10
exp_no = int(sys.argv[5])

angles_typ = sys.argv[6]
assert angles_typ in ['rand','optim'], "angles_typ should be 'rand' or 'optim'"
if angles_typ == 'rand':
    angls_typ = 'random'
if angles_typ == 'optim':
    angls_typ = 'optimized'

print(" - The model is",model)
print(" - The noise model is", noise_md_st)
print(" - ", city_no ,"cities,",layers_no,"layers,",instances_no ,"intances," "and",exp_no,"experiments")
print(" - Angles type:", angls_typ)
print()

my_keys = [None, 'yes']

ps = [[0.01, 0.005, 0.002, 0.001, 0.0005][int(sys.argv[7])]]
scheme = "0001"
# use_schem = 'on'


for p in ps:
    noise_model = nse_md(p)
    for m in range(1,exp_no+1):
        print("({}) experiment = {}".format(p, m))
        filename = "data_noise_qaoa_xy/tsp{}_plotting_for_paper/{}-model-{}-noise-{}-exp_{}-scheme{}-angles-keys_layers{}_instances-".format(city_no, model, noise_model_str,
                                                                                                                                            p, m ,scheme, layers_no)
        if os.path.exists(filename + "en.npy"):
            print("({}) > File exists - skipped")
            continue

        hamildict = np.load('./tsp_data_xy/tsp{}/tsp_dict_{}_{}.npz'.format(city_no, model, m))
        for j in [0,2]:
            for l in range(len(hamildict[:,j])):
                if hamildict[:,j][l] > 0:
                    hamildict[:,j][l] -= 1

        diag_vals = np.load('./tsp_data_xy/tsp{}/qubo_{}.npz'.format(city_no, m))
        data_en = np.zeros((len(my_keys)+1, layers_no, instances_no))
        # data_prob = np.zeros((len(my_keys)+1, layers_no, instances_no))

        for t in range(1, layers_no+1):
            print(t)
            for a in range(1,instances_no+1):
                print("({}) > {}-th instance, depth {}".format(p, a, t))
                # TODO: random angles
                if angles_typ == 'optim':
                    angle = np.load('./tsp_data_xy/data{}/qubo_{}-exp{}-m{}-result.npy'.format(city_no, model, m, a))['{}'.format(2*t)]
                    en_pure = np.load('./tsp_data_xy/data{}/qubo_{}-exp{}-m{}-result.npy'.format(city_no, model, m, a))['energies'][t-1]
                elif angles_typ == 'rand':
                    angle = np.random.random(2*t) * 2 * np.pi
                    en_pure,_ = qaoa_xy(angle,nse_md(0),diag_vals,hamiltonian=hamildict)

                #outputs energy
                # en_pure=np.load('./data_2021-04-09T11:05:44.139/data/qubo_{}-exp{}-m{}-result.npy'.format(model,m,a))['energies'][t-1]
                data_en[0,t-1,a-1] = np.real(en_pure)
                # data_prob[0,t-1,a-1] = 1
                f = lambda x: qaoa_xy(angle, noise_model,diag_vals, midmeasure = x, hamiltonian=hamildict, scheme=scheme)
                for (k_ind, k) in enumerate(my_keys):
                    data_raw = f(k)
                    data_en[k_ind+1,t-1,a-1] = np.real(data_raw[0])

        np.save(filename + "en", data_en)
        # np.save(filename + "prob", data_prob)

