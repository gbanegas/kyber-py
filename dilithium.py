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
        #if len(input_bytes) != 34:
        #    raise ValueError(f"Input bytes should be one 32 byte array and 2 single bytes.")
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
        #if len(input_bytes) != 33:
        #    raise ValueError(f"Input bytes should be one 32 byte array and one single byte.")
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
            input_bytes = self._prf(sigma,  bytes([N]), 64*self.eta)
            poly = self.R.cbd(input_bytes, self.eta, is_ntt=is_ntt)
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
            input_bytes = self._prf(sigma,  bytes([N]), 64*self.eta)
            poly = self.R.cbd(input_bytes, self.eta, is_ntt=is_ntt)
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



    def _rounding(self, t):
        print("t: {}".format(t))
        print( ((1 << (self.d-1)) -1))
        res = (t + 2048)
        print(res)
        tmm = res >> self.d
        print("tmm: {}".format(tmm))
        a1 = (t + (1 << (self.d-1) -1)) >> self.d
        a0 = t - (a1 << self.d)
        print(a1)
        return a0, a1


    def _power2Round(self, t_in):
        t1 = []
        t0 = []
        for i in range(self.k):
            tmp0 = []
            tmp1 = []
            for j in range(self.n):
                a0, a1 = self._rounding(t_in[i][0][j])
                tmp1.append(a1)
                tmp0.append(a0)
            t1.append(tmp1)
            t0.append(tmp0)

        return t1, t0

    def _polyt1_pack(self, t1):
        pk = []
        idx = int(self.n/4)
        print(t1)
        for j in range(idx):
            pk.append((t1[4*j+0] >> 0) % 255);
            pk.append(((t1[4*j+0] >> 8 ) | (t1[4*j+1] << 2 ) ) % 255);
            pk.append(((t1[4*j+1] >> 6 ) | (t1[4*j+2] << 4 )) % 255);
            pk.append(((t1[4*j+2] >> 4 ) | (t1[4*j+3] << 6 )) % 255);
            pk.append(((t1[4*j+3] >> 2 )) % 255);
        return pk


    def _pack_pk(self, rho, t1):
        pk = []
        for i in range(32):
            pk.append(rho[i])
        for i in range(self.k):
            res = self._polyt1_pack(t1[i])

            pk += res
        return pk

    def _polyeta_pack(self, pol):
        tmp = []

        if self.eta == 2:

            idx = int(self.n/8)

            for i in range(idx):
                t = []
                t.append(self.eta - pol[0][8*i+0])
                t.append(self.eta - pol[0][8*i+1])
                t.append(self.eta - pol[0][8*i+2])
                t.append(self.eta - pol[0][8*i+3])
                t.append(self.eta - pol[0][8*i+4])
                t.append(self.eta - pol[0][8*i+5])
                t.append(self.eta - pol[0][8*i+6])
                t.append(self.eta - pol[0][8*i+7])

                a = ((t[0] >> 0) | (t[1] << 3) | (t[2] << 6)) % 255
                b = ((t[2] >> 2) | (t[3] << 1) | (t[4] << 4) | (t[5] << 7)) % 255
                c = ((t[5] >> 1) | (t[6] << 2) | (t[7] << 5)) % 255

                tmp.append(a)
                tmp.append(b)
                tmp.append(c)
        if self.eta == 4:

            idx = int(self.n/2)

            for i in range(idx):
                t = []
                t.append(self.eta - pol[0][2*i+0])
                t.append(self.eta - pol[0][2*i+1])
                a = (t[0] | (t[1] << 4)) % 255
                tmp.append(a)



        return tmp

    def _polyt0_pack(self, t0):
        idx = int(self.n/8)
        tmp = []
        for i in range(idx):
            t = []
            tmp_val = (1 << (self.d-1))
            t.append(( tmp_val - t0[8*i+0]))
            t.append(( tmp_val - t0[8*i+1]))
            t.append(( tmp_val - t0[8*i+2]))
            t.append(( tmp_val - t0[8*i+3]))
            t.append(( tmp_val - t0[8*i+4]))
            t.append(( tmp_val - t0[8*i+5]))
            t.append(( tmp_val - t0[8*i+6]))
            t.append(( tmp_val - t0[8*i+7]))

            tmp.append((t[0] % 255))
        return tmp


    def _pack_sk(self, rho, key, tr, t0, s1, s2):
        sk = []
        for i in range(32):
            sk.append(rho[i])

        for i in range(32):
            sk.append(key[i])

        for i in range(48):
            sk.append(tr[i])

        for i in range(self.l):
            pack_tmp = self._polyeta_pack(s1[i])
            sk += pack_tmp

        for i in range(self.k):
            pack_tmp = self._polyeta_pack(s2[i])
            sk += pack_tmp

        for i in range(self.k):
            pack_tmp = self._polyt0_pack(t0[i])
            sk += pack_tmp



        return sk

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
        s1, N = self._generate_error_vectorl(rhoprime, N)
        s1_h = s1.to_ntt()
        s2, N = self._generate_error_vectork(rhoprime, N)
        t = A @ s1_h
        t = t.from_ntt()
        t = t + s2
        s1 = s1.from_ntt()
        t1, t0 = self._power2Round(t)
        pk = self._pack_pk(rho, t1)
        tr = list(self._kdf(bytes(pk), self.publickeybytes))
        sk = self._pack_sk(rho, key, tr, t0, s1, s2)
        #print(sk)

        return pk, sk



Dilithium2 = Dilithium(DEFAULT_PARAMETERS["dilithium_2"])

if __name__ == '__main__':
    # Test kyber_512
    pk,sk = Dilithium2.keygen()
    print("pk:\n\n {} \n\n".format(pk))
    print("sk:\n\n {} \n\n".format(sk))
