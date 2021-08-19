from encoding_qubo_decoding_mixture import *
from math import asin
from qiskit import QuantumRegister, QuantumCircuit, execute, Aer

def xy_compact(qc,theta,j,k):
    qc.cx(k, j)
    qc.crx(-1 * theta, j, k)
    qc.cx(k, j)

def convet_t1v1_to_m(t,v,n):
    m = n*(t-1)+v-1
    return m

def xy_mixer(m,qreg,theta):
    qc = QuantumCircuit(qreg)
    if m==2:
        for t in range(1,m+1):
            xy_compact(qc,theta,convet_t1v1_to_m(t,2,m),convet_t1v1_to_m(t,1,m))
    else:
        # Even cities
        for t,k in it.product(range(1,m+1),range(2,m+1,2)):
            xy_compact(qc,theta,convet_t1v1_to_m(t,k,m),convet_t1v1_to_m(t,k-1,m))
        # Odd cities
        for t,k in it.product(range(1,m+1),range(3,m+1,2)):
            xy_compact(qc,theta,convet_t1v1_to_m(t,k,m),convet_t1v1_to_m(t,k-1,m))
        # Last & 1st
        for t in range(1,m+1):
            xy_compact(qc,theta,convet_t1v1_to_m(t,1,m),convet_t1v1_to_m(t,m,m))
    # for i in range(int(m**2/2)):
    #     qc.swap(i,m**2-1-i)
    return qc

def wstate(n):
    qr = QuantumRegister(n)
    qc = QuantumCircuit(qr)
    a=n-1
    b=n
    qc = QuantumCircuit(qr)
    qc.ry(2*asin((a/b)**(1/2)),0)
    for i in range(0,n-2):
        a=n-i-2
        b=n-i-1
        qc.cry(2*asin((a/b)**(1/2)),i,i+1)
    for i in range(0,n-1):
        qc.cx(n-i-2,n-i-1)
    qc.x(0)
    return qc

def initialize_with_wstate(n,qr):
    qc = QuantumCircuit(qr)
    for m in range(1,n+1):
        qc = qc.compose(wstate(n),list(range(n*(m-1),n*(m-1)+n)))
    return qc

def classical_post_selection(n,dm):
    qw = QuantumRegister(n**2)
    w_circuit = initialize_with_wstate(n,qw)
    for i in range(int(n**2/2)):
        w_circuit.swap(i,n**2-1-i)
    backend = Aer.get_backend('statevector_simulator')
    job = execute(w_circuit, backend=backend, shots=1, memory=True)
    job_result = job.result()
    vector = job_result.get_statevector(w_circuit)
    zeros_pos = np.where(np.isclose(vector, 0, rtol=0, atol=1/np.sqrt(n**n)/2))[0]
    assert len(vector)-len(zeros_pos)==n**n, "not detecting all zero positions"
    diag_dm = np.diag(dm)
    diag_dm.setflags(write=True)
    for i in zeros_pos:
        diag_dm[i]=0
    # diag_dm = diag_dm/(sum(diag_dm))
    return diag_dm

def qaoa_xy(angles, noise_model, diag_vals, midmeasure = None, scheme=None, hamiltonian=None, normalization="full"):

    n = int(np.max(hamiltonian[:,1]))
    qr = QuantumRegister(n**2)
    main_circuit=QuantumCircuit(qr)

    times = int(len(angles)/2)
    qubo_theta = angles[0:int(times)]
    mixer_phi = angles[int(times):]

    main_circuit += initialize_with_wstate(n,qr)
    main_circuit.barrier()

    scheme = "1" if scheme == None else scheme
    for i in range(times):
        # Merging circuits
        main_circuit += qubo_circuit_simplified(hamiltonian,qr,qubo_theta[i]*2) + xy_mixer(n,qr,mixer_phi[i]*2)
        main_circuit.barrier()

        if midmeasure == 'yes' and i < times-1 and scheme[i % len(scheme)] == "1":
            main_circuit.barrier()
            main_circuit += hobo2qubo(n, qr).inverse() + kraus_gates(qr,n, True) + hobo2qubo(n, qr)
            main_circuit.barrier()

    # Swap all
    for i in range(int(n**2/2)):
        main_circuit.swap(i,n**2-1-i)

    # Measure Energy
    backend = QasmSimulator(method='density_matrix', noise_model=noise_model,
                            basis_gates = ['cx', 'u1','u2','u3','kraus'])
    main_circuit.append(Snapshot('ss', 'density_matrix', n**2), qr)
    job = execute(main_circuit, backend=backend)
    result = job.result().results[0].to_dict()['data']['snapshots']['density_matrix']['ss'][0]['value']
    dm = np.asmatrix(result)

    diag_dm = classical_post_selection(n,dm)
    prob = np.sum(diag_dm)

    #find mean energy
    mean_energy = np.real(np.sum(diag_dm * diag_vals))/prob
    max_energy = max(diag_vals)
    min_energy = min(diag_vals)
    # energy = mean_energy
    if normalization == "full":
        energy = (mean_energy-min_energy)/(max_energy-min_energy)
    elif normalization == "shift":
        energy = mean_energy-min_energy

    # return energy,prob
    return energy,prob
