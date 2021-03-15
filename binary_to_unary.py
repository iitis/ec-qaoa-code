from qiskit import QuantumRegister, QuantumCircuit, ClassicalRegister, execute, BasicAer, Aer
import numpy as np
import math as math

def binary_to_unary(qreg, n, binary_bits=None, decomposition = 'on'): # TODO binary_bits should have some default input
    assert n > 1
    kmax = int(math.ceil(math.log(n, 2)))
    if binary_bits == None:
        binary_bits = [2**k-1 for k in range(1, kmax)]
    binary_bits.append(n-1)
    assert qreg.size == n # assert number of qubits in qreg is exactly N
    qc = QuantumCircuit(qreg)

    qc.x(0)
    qc.cx(1, 0)
    if n == 2:
        return qc, binary_bits

    for logk in range(2, kmax):
        k = 2**logk
        kprev = int(k/2)
        for i in range(0, kprev-1):
            if decomposition == 'on':
                qc.cx(i+kprev,i)
                qc.rccx(n-1,i,i+kprev)
                qc.cx(i+kprev,i)
            elif decomposition == 'off':
                qc.cswap(k-1, i, i + kprev)

        for i in range(kprev, k-1):
            qc.cx(i, k-1)
        qc.cx(k-1, kprev -1)

    kprev = 2**(kmax-1)
    for i in range(0, n-1-kprev):
        if decomposition == 'on':
            qc.cx(i+kprev,i)
            qc.rccx(n-1,i,i+kprev)
            qc.cx(i+kprev,i)
        elif decomposition == 'off':
            qc.cswap(n-1, i, i + kprev)

    for i in range(kprev, n-1):
        qc.cx(i, n-1)
    qc.cx(n-1, kprev -1 - (2**kmax-n))
    return qc, binary_bits

if __name__ == "__main__":
    for i in range(3,4):
        hamildict = np.load('./qaoa_efficiency_analysis/tsp_dict/tsp_{}-1.npz'.format(i))
        n = int(np.max(hamildict[:,1]))
        qr = QuantumRegister(n)
        qc, binary_bits = binary_to_unary(qr, n)
        depth = qc.depth()
        print(n+1,'->',depth)
