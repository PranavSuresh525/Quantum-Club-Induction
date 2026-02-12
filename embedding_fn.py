from qiskit import QuantumCircuit
import numpy as np
def create_embedding_circuit(n, secret_image, intensity_bits=8):
    size=2**n
    total_qubits=2*n+intensity_bits

    coord_qubits=list(range(0, 2*n))
    intensity_qubits=list(range(2*n, total_qubits))
    lsb_qubit= total_qubit-1

    qc=QuantumCircuit(total_qubits)

    def embedding_fn( x_val, y_val):
        x_bit=format(x_val, f'0{n}b')
        y_bit=format(y_val, f'0{n}b')

        full_bits= x_bits+ y_bits

        for i, bit in enumerate(full_bits):
            if bit==0:
                 qc.x(coord_qubit[i])
        qc.mcx(coord_bit, lsb_qubit)

        for i, bit in enumerate( full_bits):
            if bit==0:
               qc.x(coord_qubit[i])

    for x in range(size):
        for y in range(size):
            if secret_image[x,y]==1:
                embedding_fn(x,y)

    return qc               
