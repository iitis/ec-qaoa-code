from qiskit import QuantumRegister,ClassicalRegister,QuantumCircuit,assemble,Aer,execute, IBMQ
import numpy as np
import math
import itertools as it

def convet_t1v1_to_m(t,v,n):
    m = n*(t-1)+v-1
    return m

def qubo_circuit_simplified(hamildict,q,theta):
    # Initialize quantum circuit
    qc = QuantumCircuit(q)
    n = int(np.max(hamildict[:,1]))

    int_dict={tuple(map(int, x[0:4])) : x[4] * theta for x in hamildict}

    def double_gates(m1,m2,w):
        qc.cx(q[m1],q[m2])
        qc.rz(w,q[m2])
        qc.cx(q[m1],q[m2])

    for i in int_dict.keys():
        if (i[0]!=0) and (i[2]!=0):
            t1 = i[0]
            v1 = i[1]
            t2 = i[2]
            v2 = i[3]
            w = int_dict[i]
            m1 = convet_t1v1_to_m(t1,v1,n)
            m2 = convet_t1v1_to_m(t2,v2,n)
            double_gates(m1,m2,w)
        if i[0] != 0 and i[2] == 0:
            t = i[0]
            v = i[1]
            w = int_dict[(t,v,0,0)]
            qc.rz(w,q[convet_t1v1_to_m(t,v,n)])
    return qc

if __name__ == "__main__":
    theta = 0
    hamildict = np.load('./qaoa_efficiency_analysis/tsp_dict/tsp_5-1.npz')
    n = int(np.max(hamildict[:,1]))
    qr = QuantumRegister(n**2,'q')
    qc = qubo_circuit_simplified(hamildict,qr,theta)
    depth = qc.depth()

    print(depth/n)
    # print(qc)
