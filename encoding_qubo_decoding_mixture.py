from binary_to_unary import binary_to_unary
from cpkraus import CPKraus
from matplotlib import rc
from qiskit import *
from qiskit.providers.aer.extensions import Snapshot
from qiskit.providers.aer import QasmSimulator
import qiskit.providers.aer.noise as noise
from qiskit.quantum_info import Kraus,Statevector, DensityMatrix
from qiskit.quantum_info.operators.base_operator import BaseOperator
from qiskit.quantum_info.operators import Operator, Pauli
from qiskit.tools.visualization import plot_histogram, plot_bloch_multivector
from qubo_circuit_sparse import *
from roundRobin import roundRobin
import itertools as it
import matplotlib.pyplot as plt
import numpy as np
# %matplotlib inline #TODO: work on the visualization

# Unary to binary implementation
def hobo2qubo(n, qr):
    qreg = QuantumRegister(n)
    qc_ub, _ = binary_to_unary(qreg, n)
    qc_init = QuantumCircuit(qr)
    init_list = list(range(n))
    new_qc = qc_init + qc_init.compose(qc_ub,init_list)
    for i in range (n-1):
        init_list = [list + n for list in init_list]
        new_qc = new_qc + qc_init.compose(qc_ub,init_list)
    return new_qc

# Position of Hadamard Gates:
def hadamards_position(n,qr):
    qreg = QuantumRegister(n)
    qc = QuantumCircuit(qr)
    _, binary_bits = binary_to_unary(qreg, n)
    for i in range(1,n+1):
        for k in binary_bits:
            qc.h(n*(i-1)+k)
    return qc

# Mixer
def make_mixer(n,qr,phi):
    qreg = QuantumRegister(n)
    qc = QuantumCircuit(qr)
    _, binary_bits = binary_to_unary(qreg, n)
    for i in range(1,n+1):
        for k in binary_bits:
            qc.rx(phi,n*(i-1)+k)
    return qc

# Kraus Gates
def kraus_gates(qr,n, postselection = True):
    qc = QuantumCircuit(qr)
    if postselection:
        chan = CPKraus([[[1, 0],[ 0, 0]], [[0, 0],[ 0, 0]]])
    else:
        chan = CPKraus([[[1, 0],[ 0, 0]], [[0, 0],[ 0, 1]]])
    res_list = []
    l = [0,0]
    qreg = QuantumRegister(n)
    _, binary_bits = binary_to_unary(qreg, n)
    for i in range(n):
        zipped_lists = zip(binary_bits, l)
        s = [x + y for (x, y) in zipped_lists]
        l1 = zip(l, [n,n])
        l = [x + y for (x, y) in l1]
        res_list += s
    for q in range(0, n**2):
        if q not in res_list:
            qc.append(chan, [q])
    return qc

# Measurement part
def hobo_qubo_mixer_measurement(hamildict, noise_model,angles,diag_vals, ps_end=False, mid_measure = None):

    n = int(np.max(hamildict[:,1]))
    qr = QuantumRegister(n**2)
    main_circuit=QuantumCircuit(qr)

    times = len(angles)
    qubo_theta = angles[0:int(times/2)]
    mixer_phi = angles[int(times/2):]
    main_circuit += hadamards_position(n,qr)

    # Merging circuits
    for theta,phi in zip(qubo_theta,mixer_phi):
        main_circuit += hobo2qubo(n, qr) + qubo_circuit_simplified(hamildict,qr,theta) + hobo2qubo(n, qr).inverse()
        main_circuit += make_mixer(n,qr,phi)
        
        if mid_measure == 'no-ps':
            main_circuit += kraus_gates(qr,n, False)
        elif mid_measure == 'ps':
            main_circuit += kraus_gates(qr,n, True)
        else:
            pass

    if ps_end == 'end':
        main_circuit += kraus_gates(qr,n, True)

    main_circuit += hobo2qubo(n, qr)
    
    # Swap all
    for i in range(int(n**2/2)):
        main_circuit.swap(i,n**2-1-i)

    #Preparation of Density-Matrix:
    backend = QasmSimulator(method='density_matrix', noise_model=noise_model,
                            basis_gates = ['cx', 'u1','u2','u3','kraus'])
    main_circuit.append(Snapshot('ss', 'density_matrix', n**2), qr)
    job = execute(main_circuit, backend=backend)
    result = job.result().results[0].to_dict()['data']['snapshots']['density_matrix']['ss'][0]['value']
    dm = np.asmatrix(result)

    #finds probability
    prob = np.trace(dm.real)

    #find mean energy
    mean_energy = np.real(np.sum(np.diag(dm) * diag_vals))/prob
    max_energy = max(diag_vals)
    min_energy = min(diag_vals)
    energy = (mean_energy-min_energy)/(max_energy-min_energy)
    return energy, prob

def hobo_qubo_mixer_measurement_pure(hamildict,angles,diag_vals):

    n = int(np.max(hamildict[:,1]))
    qr = QuantumRegister(n**2)
    main_circuit=QuantumCircuit(qr)

    times = len(angles)
    qubo_theta = angles[0:int(times/2)]
    mixer_phi = angles[int(times/2):]
    main_circuit += hadamards_position(n,qr)

    # Merging circuits
    for theta,phi in zip(qubo_theta,mixer_phi):
        main_circuit += hobo2qubo(n, qr) + qubo_circuit_simplified(hamildict,qr,theta) + hobo2qubo(n, qr).inverse()
        main_circuit += make_mixer(n,qr,phi)

    main_circuit += hobo2qubo(n, qr)
    
    # Swap all
    for i in range(int(n**2/2)):
        main_circuit.swap(i,n**2-1-i)

    #Preparation of Density-Matrix:
    backend = Aer.get_backend('statevector_simulator')
    result = execute(main_circuit, backend).result()
    statevector = result.get_statevector(main_circuit)

    #find mean energy
    mean_energy = np.real(np.sum(np.abs(statevector)**2 * diag_vals))
    max_energy = max(diag_vals)
    min_energy = min(diag_vals)
    energy = (mean_energy-min_energy)/(max_energy-min_energy)
    return energy, 1.


if __name__ == "__main__":
    hamildict = np.load('./qaoa_efficiency_analysis/tsp_dict/tsp_3-1.npz')
    angles = np.load('./qaoa_efficiency_analysis/tsp_dict/angles_test_3-1.npz')
    diag_vals = np.load('./qaoa_efficiency_analysis/tsp_dict/qubo_vec_3-1.npz')
    noise_model = create_noise_model(0,0)
    mean_energy = hobo_qubo_mixer_measurement(hamildict,noise_model,angles,diag_vals,'end')
    print(main_circuit)
    print('Hooray! the mean energy is', mean_energy)