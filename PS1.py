import numpy as np

class Quantum_register:
    X = np.array([[0,1],[1,0]])
    H = np.array([[1,1],[1,-1]])/np.sqrt(2)
    Y = np.array([[0,-1j],[1j,0]])
    Z = np.array([[1,0],[0,-1]])
    S = np.array([[1,0],[0,1j]])
    T = np.array([[1,0],[0,np.exp(1j*np.pi/4)]])
    CX = np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
    CZ = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,-1]])
    CY = np.array([[1,0,0,0],[0,1,0,0],[0,0,0,-1j],[0,0,1j,0]])
    CNOT = CX
    SWAP = np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])

    def __init__(self, n_qubits):
        self.n = n_qubits
        self.statevector = np.zeros((2**n_qubits,1),dtype=complex)
        self.statevector[0][0] = 1.0

    def apply_gate(self, gate, qubits):
        if len(qubits) == 1:
            I = np.eye(2)
            U = np.kron(gate, I) if qubits[0]==0 else np.kron(I, gate)
            self.statevector = U @ self.statevector
        elif len(qubits) == 2:
            ctrl, target = qubits
            if [ctrl, target] == [0,1]:
                self.statevector = gate @ self.statevector
            elif [ctrl, target] == [1,0]:
                SWAP = Quantum_register.SWAP
                self.statevector = SWAP @ gate @ SWAP @ self.statevector

    def measure(self):
        probs = np.abs(self.statevector)**2
        probs /= np.sum(probs)
        outcome = np.random.choice(2**self.n, p=probs.flatten())
        return format(outcome, f'0{self.n}b')

qr = Quantum_register(2)
qr.apply_gate(Quantum_register.X, [1])
qr.apply_gate(Quantum_register.CX, [0,1])
print(qr.statevector)
print(qr.measure())
