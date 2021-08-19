import qiskit.providers.aer.noise as noise
import numpy as np

# amplitude damping
def create_noise_model_amp_damping(p):
    error = noise.amplitude_damping_error(p)
    error2 = error.tensor(error)
    noise_model = noise.NoiseModel()
    noise_model.add_all_qubit_quantum_error(error, ['u1', 'u2', 'u3'])
    noise_model.add_all_qubit_quantum_error(error2, ['cx'])
    return noise_model

# depolarizing channel
def create_noise_model_depolarazing(p1,p2=None):
    error_1 = noise.depolarizing_error(p1, 1)
    if p2==None:
        p2 = p1
    error_2 = noise.depolarizing_error(p2, 2)
    noise_model = noise.NoiseModel()
    noise_model.add_all_qubit_quantum_error(error_1, ['u1', 'u2', 'u3'])
    noise_model.add_all_qubit_quantum_error(error_2, ['cx'])
    return noise_model

# random unitary noise
def create_noise_model_random_x(p):
    u = [[0,  1], [1, 0]]
    error1 = noise.mixed_unitary_error([(u, p), (np.eye(2),1-p)])
    error2 = error1.tensor(error1)
    noise_model = noise.NoiseModel()
    noise_model.add_all_qubit_quantum_error(error1, ['u1', 'u2', 'u3'])
    noise_model.add_all_qubit_quantum_error(error2, ['cx'])
    return noise_model

if __name__ == "__main__":
    noise_model = create_noise_model_amp_damping(0.5)
