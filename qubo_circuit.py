# Testing circuit from arxiv:2009.07309
from qiskit import QuantumRegister,ClassicalRegister,QuantumCircuit,assemble,Aer,execute, IBMQ
import numpy as np
import math
import itertools as it

def convet_t1v1_to_m(t,v,n):
    m = n*(t-1)+v-1
    return m

def roundRobin(units, sets=None):
    """ Generates a schedule of "fair" pairings from a list of units """
    if len(units) % 2:
        units.append(None)
    count    = len(units)
    sets     = sets or (count - 1)
    half     = count / 2
    schedule = []
    for turn in range(sets):
        pairings = []
        for i in range(int(half)):
            pairings.append(units[i])
            pairings.append(units[count-i-1])
        units.insert(1, units.pop())
        schedule.append(pairings)
    return schedule

def check_none(b):
    return None not in b

def roundrobinlist(n):
    cities = list(range(1,n+1))
    tournament = []
    for pairings in roundRobin(cities):
        tournament.append(list(zip(pairings[0::2],pairings[1::2])))
    newlist=[]
    for x in tournament:
        newlist.append(list(filter(check_none,x)))
    return newlist

def qubo_circuit (hamildict,q,theta):
    # Initialize quantum circuit
    qc = QuantumCircuit(q)
    n = int(np.max(hamildict[:,1]))

    int_dict={tuple(map(int, x[0:4])) : x[4] * theta for x in hamildict}

    def double_gates(m1,m2,w):
        qc.cx(q[m1],q[m2])
        qc.rz(w,q[m2])
        qc.cx(q[m1],q[m2])

    # Round Robin Blocks
    newlist=roundrobinlist(n)
    def rrt_blocks(block_typ):
        for tv, r in it.product(list(range(1,n+1)),newlist):
            for match_x in r:
                match = sorted(match_x)
                if block_typ == 1:
                    t1 = t2 = tv
                    v1 = match[0]
                    v2 = match[1]
                elif block_typ == 2:
                    t1 = match[0]
                    t2 = match[1]
                    v1 = v2 = tv
                else:
                    print(' Please incert 1, for block A or 2 for block B')
                w = int_dict[(t1,v1,t2,v2)]
                m1 = convet_t1v1_to_m(t1,v1,n)
                m2 = convet_t1v1_to_m(t2,v2,n)
                double_gates(m1,m2,w)

    # - Block A
    rrt_blocks(1)
    qc.barrier()

    # - Block B

    rrt_blocks(2)
    qc.barrier()

    # C Blocks
    def c_block(block_typ):
        for t, k, i in it.product(range(block_typ,n,2),range(1,n),range(1,n)):
            t1 = t
            t2 = t+1
            if t2 > n:
                t2 -= (n)
            v1 = i
            v2 = i+k
            if v2 > n:
                v2 -= (n)
            if (t1, v1, t2, v2) not in int_dict.keys():
                t1,v1,t2,v2 = t2,v2,t1,v1

            w = int_dict[(t1, v1, t2, v2)]
            m1 = convet_t1v1_to_m(t1,v1,n)
            m2 = convet_t1v1_to_m(t2,v2,n)
            double_gates(m1,m2,w)

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