from qiskit import QuantumRegister, QuantumCircuit, ClassicalRegister, execute, BasicAer, Aer
import numpy as np
import math as math

def binary_to_unary(qreg, n, binary_bits=None): 
    assert n > 1
    kmax = int(math.ceil(math.log(n, 2)))
    if binary_bits == None:
        binary_bits = [2**k-1 for k in range(1, kmax)]
    binary_bits.append(n-1)
    assert qreg.size == n
    qc = QuantumCircuit(qreg) 

    qc.x(0)
    qc.cx(1, 0)
    qc.barrier()
    if n == 2:
        return qc, binary_bits

    for logk in range(2, kmax):
        k = 2**logk
        kprev = int(k/2)
        for i in range(0, kprev-1):
            qc.cswap(k-1, i, i + kprev)
        for i in range(kprev, k-1):
            qc.cx(i, k-1)
        qc.cx(k-1, kprev -1)
        qc.barrier()

    kprev = 2**(kmax-1)
    for i in range(0, n-1-kprev):
        qc.cswap(n-1, i, i + kprev)
    for i in range(kprev, n-1):
        qc.cx(i, n-1)
    qc.cx(n-1, kprev -1 - (2**kmax-n))
    return qc, binary_bits


def relevant_positions_binary(n, binary_bits):
    qreg = QuantumRegister(n)
    qc = QuantumCircuit(qreg)
    for b in binary_bits:
        qc.h(b)
    backend = BasicAer.get_backend('statevector_simulator')   
    job = execute(qc, Aer.get_backend('statevector_simulator'))
    the_state = job.result().get_statevector(qc)
    positions = list(filter(lambda i:np.abs(the_state[i])**2>1/(2*n), range(len(the_state))))
    
    assert len(positions) == 2**int(math.ceil(math.log(n, 2)))
    return positions

def test(n):
    backend = BasicAer.get_backend('qasm_simulator')
    qreg = QuantumRegister(n)
    creg = ClassicalRegister(n)
    qc, binary_bits = binary_to_unary(qreg, n)
    print(binary_bits)
    for i in range(n):
        qc_test = QuantumCircuit(qreg, creg)
        for bits, pos in zip(reversed(bin(i)), binary_bits): 
            if bits == '1':
                qc_test.x(pos)
        qc_test.barrier()
        qc_test += qc 
        qc_test.measure(qreg, creg)

        result_sim = execute(qc_test, backend, shots=1024).result().get_counts(qc_test)
        keys = result_sim.keys()
        print(keys)

if __name__ == "__main__":
    for n in range(2,4):
        qreg = QuantumRegister(n)
        qc, binary_bits = binary_to_unary(qreg, n)
        pos = relevant_positions_binary(n, binary_bits)
        
        qc1 = QuantumCircuit(qreg)
        for b in binary_bits:
            qc1.h(b)
        qc1.barrier()
        print((qc1+qc).draw())

        job = execute(qc, Aer.get_backend('unitary_simulator'))
        u = np.real(job.result().get_unitary(qc)[:,pos])
        np.save("permutations_cutters/permcut_{}".format(n), u.transpose())
