import os
from hashlib import sha3_256, sha3_512, shake_128, shake_256
from polynomials import *
from modules import *

DEFAULT_PARAMETERS = {
    "dilithium_2" : {
        "n" : 256,
        "k" : 4,
        "q" : 8380417,
        "root_unity": 1753,
        "d": 13,
        "eta" : 2,
        "l" : 4,
        "tau" : 39,
        "beta": 78,
        "omega": 80,
        "gamma1": 131072,
        "gamma2": 95232,
        "publickeybytes": 1312,
    },
    "dilithium_3" : {
        "n" : 256,
        "k" : 3,
        "q" : 8380417,
        "eta_1" : 2,
        "eta_2" : 2,
        "du" : 10,
        "dv" : 4,
        "root_unity": 1753,
        "d": 13,
    },
    "dilithium_5" : {
        "n" : 256,
        "k" : 4,
        "q" : 8380417,
        "eta_1" : 2,
        "eta_2" : 2,
        "du" : 11,
        "dv" : 5,
        "root_unity": 1753,
        "d": 13,
    }
}

class Dilithium:
    def __init__(self, parameter_set):
        self.n = parameter_set["n"]
        self.k = parameter_set["k"]
        self.q = parameter_set["q"]
        self.eta = parameter_set["eta"]
        self.root_unity = parameter_set["root_unity"]
        self.d = parameter_set["d"]
        self.l = parameter_set["l"]
        self.tau = parameter_set["tau"]
        self.beta = parameter_set["beta"]
        self.omega = parameter_set["omega"]
        self.gamma1 = parameter_set["gamma1"]
        self.gamma2 = parameter_set["gamma2"]
        self.publickeybytes = parameter_set["publickeybytes"]

        self.R = PolynomialRing(self.q, self.n)
        self.M = Module(self.R)

    @staticmethod
    def _xof(bytes32, a, b, length):
        """
        XOF: B^* x B x B -> B*
        """
        input_bytes = bytes32 + a + b
        if len(input_bytes) != 34:
            raise ValueError(f"Input bytes should be one 32 byte array and 2 single bytes.")
        return shake_128(input_bytes).digest(length)

    @staticmethod
    def _h(input_bytes):
        """
        H: B* -> B^32
        """
        return sha3_256(input_bytes).digest()

    @staticmethod
    def _g(input_bytes):
        """
        G: B* -> B^32 x B^32
        """
        output = sha3_512(input_bytes).digest()
        return output[:32], output[32:]

    @staticmethod
    def _prf(s, b, length):
        """
        PRF: B^32 x B -> B^*
        """
        input_bytes = s + b
        if len(input_bytes) != 33:
            raise ValueError(f"Input bytes should be one 32 byte array and one single byte.")
        return shake_256(input_bytes).digest(length)

    @staticmethod
    def _kdf(input_bytes, length):
        """
        KDF: B^* -> B^*
        """
        return shake_256(input_bytes).digest(length)

    def _generate_error_vectorl(self, sigma, N, is_ntt=False):
        """
        Helper function which generates a element in the
        module from the Centered Binomial Distribution.
        """
        elements = []
        for i in range(self.l):
            input_bytes = self._prf(sigma,  bytes([N]), 64*self.eta_1)
            poly = self.R.cbd(input_bytes, self.eta_1, is_ntt=is_ntt)
            elements.append(poly)
            N = N + 1
        v = self.M(elements).transpose()
        return v, N

    def _generate_error_vectork(self, sigma, N, is_ntt=False):
        """
        Helper function which generates a element in the
        module from the Centered Binomial Distribution.
        """
        elements = []
        for i in range(self.k):
            input_bytes = self._prf(sigma,  bytes([N]), 64*self.eta_1)
            poly = self.R.cbd(input_bytes, self.eta_1, is_ntt=is_ntt)
            elements.append(poly)
            N = N + 1
        v = self.M(elements).transpose()
        return v, N

    def _generate_matrix_from_seed(self, rho, N, transpose=False, is_ntt=False):
        """
        Helper function which generates a element of size
        k x l from a seed `rho`.

        When `transpose` is set to True, the matrix A is
        built as the transpose.
        """
        A = []
        for i in range(self.k):
            row = []
            for j in range(self.l):
                if transpose:
                    input_bytes = self._xof(rho, bytes([i]), bytes([j]), 3*self.R.n)
                else:
                    input_bytes = self._xof(rho, bytes([j]), bytes([i]), 3*self.R.n)
                aij = self.R.parse(input_bytes, is_ntt=is_ntt)
                row.append(aij)
                N = N + 1
            A.append(row)
        return self.M(A), N


    @staticmethod
    def rounding(self, t):
        a1 = (a+ (1 << (self.d-1) -1)) >> self.d
        a0 = a - (a1 << self.d)
        return a0, a1

    @staticmethod
    def _power2Round(self, t):
        t1 = []
        t0 = []
        for i in range(self.k):
            for j in range(self.n):
                a0, a1 = rounding(t[i][j])
                t1.append(a1)
                t0.append(a0)
        return t1, t0

    @staticmethod
    def pack_pk(self, rho, t1):
        t1enc = t1.encode()
        return rho+t1enc

    def keygen(self):
        """


        Output:
            pk: Public key
            sk: Secret key

        """
        r1 = os.urandom(32)
        rho = self._kdf(r1, (2*32)+64)
        rhoprime = rho[32::]
        key = rhoprime[64::]
        N = 0
        A, N = self._generate_matrix_from_seed(rho, N, transpose=False, is_ntt=True)
        s1, N = self._generate_error_vectorl(rhoprime, N).to_ntt()
        s2, N = self._generate_error_vectork(rhoprime, N)
        t = A @ s1
        t = t.from_ntt()
        t = t + s2
        t1, t0 = self._power2Round(t)
        pk = delf._pack_pk(rho, t1)
        tr = self._kdf(pk, self.publickeybytes)

        return pk



Dilithium2 = Dilithium(DEFAULT_PARAMETERS["dilithium_2"])

if __name__ == '__main__':
    # Test kyber_512
    pk = Dilithium2.keygen()
