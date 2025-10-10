import cmath, math, random

class np:
    pi = math.pi
    @staticmethod
    def sqrt(x): return math.sqrt(x)
    @staticmethod
    def exp(x): return cmath.exp(x)
    @staticmethod
    def eye(n):
        m = [0j] * (n * n)
        for i in range(n):
            m[i * n + i] = 1.0
        return m
    @staticmethod
    def kron(A, B):
        nA = int(math.sqrt(len(A)))
        nB = int(math.sqrt(len(B)))
        out = [0j] * (nA * nB * nA * nB)
        for i in range(nA):
            for j in range(nA):
                for k in range(nB):
                    for l in range(nB):
                        out[(i * nB + k) * (nA * nB) + (j * nB + l)] = A[i * nA + j] * B[k * nB + l]
        return out
    class random:
        @staticmethod
        def choice(n, p):
            r = random.random()
            total = 0
            for i, prob in enumerate(p):
                total += prob
                if r < total:
                    return i
            return n - 1

class Quantum_register:
    def __init__(self, n_qubits):
        self.n = n_qubits
        self.statevector = [0j] * (2 ** n_qubits)
        self.statevector[0] = 1.0
    @staticmethod
    def matmul(A, v):
        n = int(math.sqrt(len(A)))
        res = [0j] * n
        for i in range(n):
            for j in range(n):
                res[i] += A[i * n + j] * v[j]
        return res
    @staticmethod
    def matmul_square(A, B):
        n = int(math.sqrt(len(A)))
        res = [0j] * (n * n)
        for i in range(n):
            for j in range(n):
                res[i * n + j] = sum(A[i * n + k] * B[k * n + j] for k in range(n))
        return res
    X = [0, 1, 1, 0]
    H = [1 / math.sqrt(2), 1 / math.sqrt(2), 1 / math.sqrt(2), -1 / math.sqrt(2)]
    Y = [0, -1j, 1j, 0]
    Z = [1, 0, 0, -1]
    S = [1, 0, 0, 1j]
    T = [1, 0, 0, cmath.exp(1j * math.pi / 4)]
    CX = [1,0,0,0, 0,1,0,0, 0,0,0,1, 0,0,1,0]
    CZ = [1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,-1]
    CY = [1,0,0,0, 0,1,0,0, 0,0,0,-1j, 0,0,1j,0]
    CNOT = CX
    SWAP = [1,0,0,0, 0,0,1,0, 0,1,0,0, 0,0,0,1]
    def apply_gate(self, gate, qubits):
        if len(qubits) == 1:
            if self.n == 1:
                self.statevector = self.matmul(gate, self.statevector)
            else:
                I = np.eye(2)
                U = np.kron(gate, I) if qubits[0] == 0 else np.kron(I, gate)
                self.statevector = self.matmul(U, self.statevector)
        elif len(qubits) == 2:
            ctrl, target = qubits
            if [ctrl, target] == [0, 1]:
                self.statevector = self.matmul(gate, self.statevector)
            elif [ctrl, target] == [1, 0]:
                SWAP = Quantum_register.SWAP
                U = self.matmul_square(SWAP, gate)
                U = self.matmul_square(U, SWAP)
                self.statevector = self.matmul(U, self.statevector)
    def measure(self):
        probs = [abs(a) ** 2 for a in self.statevector]
        total = sum(probs)
        probs = [p / total for p in probs]
        outcome = np.random.choice(2 ** self.n, p=probs)
        return format(outcome, f'0{self.n}b')

print("---- 1-Qubit Example ----")
qr1 = Quantum_register(1)
qr1.apply_gate(Quantum_register.H, [0])
print("Statevector:", qr1.statevector)
print("Measurement:", qr1.measure())

print("\n---- 2-Qubit Example ----")
qr2 = Quantum_register(2)
qr2.apply_gate(Quantum_register.X, [1])
qr2.apply_gate(Quantum_register.CX, [0, 1])
print("Statevector:", qr2.statevector)
print("Measurement:", qr2.measure())
