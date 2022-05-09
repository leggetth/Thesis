import math

import numpy
import scipy.special
import scipy.fftpack
import nistmatrixdefinitions as nmatrix
from random import choice
import scipy.constants

def freq_monobit(bitstring):
    array = [int(x) for x in bitstring]
    ones = array.count(1)
    zeros = array.count(0)
    difference = abs(ones - zeros)
    s_obs = float(difference) / math.sqrt(float(len(bitstring)))
    p_value = math.erfc(s_obs / math.sqrt(2))
    return p_value

def block_freq(block_len, bitstring):
    num_block = len(bitstring) // block_len
    sum = 0
    for i in range(num_block):
        block_sum = 0
        for j in range(num_block):
            block_sum += int(bitstring[j+i*num_block])
        pi = block_sum / num_block
        v = pi - 0.5
        sum += v * v
    chi_squared = 4.0 * block_len * sum
    p_value = scipy.special.gammaincc((num_block/2.0), (chi_squared/2.0))
    return p_value

def runs(bitstring):
    n = len(bitstring)
    array = [int(x) for x in bitstring]
    S = array.count(1)
    pi = float(S) / n
    if abs(pi - 0.5) > (2.0 / math.sqrt(n)):
        p_value = 0.0
    else:
        V = 1
        for k in range(1, n):
            if array[k] != array[k-1]:
                V += 1
        erfc_arg = abs(V - 2.0 * n * pi * (1 - pi)) / (2.0 * pi * (1 - pi) * math.sqrt(2 * n))
        p_value = math.erfc(erfc_arg)
    return p_value

def longestrunofones(bitstring):
    n = len(bitstring)
    assert n >= 128
    if n < 6272:
        K = 3
        M = 8
        V = [1, 2, 3, 4]
        pi = [0.21484375, 0.3671875, 0.23046875, 0.1875]
    elif 6272 < n < 750000:
        K = 5
        M = 128
        V = [4, 5, 6, 7, 8, 9]
        pi = [0.1174035788, 0.242955959, 0.249363483, 0.17517706, 0.102701071, 0.112398847]
    else:
        K = 6
        M = 10000
        V = [10, 11, 12, 13, 14, 15, 16]
        pi = [0.0882, 0.2092, 0.2483, 0.1933, 0.1208, 0.0675, 0.0727]

    N = n // M
    nu = [0] * (K + 1)
    for i in range(N):
        v_n_obs = 0
        run = 0
        for j in range(M):
            if int(bitstring[i*M+j]) == 1:
                run += 1
                if run > v_n_obs:
                    v_n_obs = run
            else:
                run = 0
        if v_n_obs < V[0]:
            nu[0] += 1
        for j in range(K+1):
            if v_n_obs == V[j]:
                nu[j] += 1
        if v_n_obs > V[K]:
            nu[K] += 1

    chi2 = 0.0
    for i in range(K+1):
        chi2 += ((nu[i] - N * pi[i]) * (nu[i] - N * pi[i])) / (N * pi[i])

    p_value = scipy.special.gammaincc(float(K/2.0), (chi2 / 2.0))
    return p_value

def binary_matrix_rank(bitstring):
    n = len(bitstring)
    M = int(math.sqrt(n / 38))
    Q = int(math.sqrt(n / 38))
    N = n // (M * Q)

    if N == 0:
        p_value = 0.0
    else:
        r = M
        product = 1
        for i in range(r):
            product *= float(((1.0 - (2 ** (i - Q))) * (1.0 - (2 ** (i - M)))) / (1.0 - (2 ** (i - r))))
        p_32 = float(2 ** (r * (Q + M - r) - (M * Q))) * product


        r = M - 1
        product = 1
        for i in range(r):
            product *= float(((1.0 - (2 ** (i - Q))) * (1.0 - (2 ** (i - M)))) / (1.0 - (2 ** (i - r))))
        p_31 = float(2 ** (r * (Q + M - r) - (M * Q))) * product

        print(p_31, p_32)
        p_30 = 1 - (p_32 + p_31)

        F_32 = 0
        F_31 = 0
        for k in range(N):
            matrix = nmatrix.BinaryMatrix(bitstring, M, Q, k)
            R = matrix.compute_rank()
            if R == M:
                F_32 += 1
            if R == M - 1:
                F_31 += 1
        F_30 = (N - (F_32 + F_31))
        chi_squared_F_32 = ((F_32 - (N * p_32)) ** 2) / (N * p_32)
        chi_squared_F_31 = ((F_31 - (N * p_31)) ** 2) / (N * p_31)
        chi_squared_F_30 = ((F_30 - (N * p_30)) ** 2) / (N * p_30)
        chi_squared = chi_squared_F_32 + chi_squared_F_31 + chi_squared_F_30
        arg1 = -chi_squared/2.0
        p_value = math.e ** arg1
    return p_value

def discrete_fourier_transform(bitstring):
    n = len(bitstring)
    X = [0] * n
    m = [0] * (n // 2)
    for i in range(n):
        X[i] = 2 * int(bitstring[i]) - 1

    dft = scipy.fftpack.rfft(X, n=n)

    m[0] = math.sqrt(dft[0] * dft[0])

    i = 0
    while i < (n // 2) - 1:
        m[i + 1] = math.sqrt((dft[2*i+1] ** 2) + (dft[2*i+2] ** 2))
        i += 1

    count = 0
    upperbound = math.sqrt(math.log(1.0/0.05) * n)
    for i in range(n // 2):
        if m[i] < upperbound:
            count += 1

    percentile = count / (n / 2) * 100
    N_1 = count
    N_0 = 0.95 * n / 2.0
    d = (N_1 - N_0) / math.sqrt(n / 4.0 * 0.95 * 0.05)
    p_value = math.erfc(abs(d) / math.sqrt(2.0))
    return p_value


def nonOverlappingTemplateMatchings(bitstring):
    n = len(bitstring)
    maxnumberOfTemplates = 148
    K = 5
    m = 9

    templates_m9 = [[0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 1, 1], [0, 0, 0, 0, 0, 0, 1, 0, 1],
                    [0, 0, 0, 0, 0, 0, 1, 1, 1], [0, 0, 0, 0, 0, 1, 0, 0, 1], [0, 0, 0, 0, 0, 1, 0, 1, 1],
                    [0, 0, 0, 0, 0, 1, 1, 0, 1], [0, 0, 0, 0, 0, 1, 1, 1, 1], [0, 0, 0, 0, 1, 0, 0, 0, 1],
                    [0, 0, 0, 0, 1, 0, 0, 1, 1], [0, 0, 0, 0, 1, 0, 1, 0, 1], [0, 0, 0, 0, 1, 0, 1, 1, 1],
                    [0, 0, 0, 0, 1, 1, 0, 0, 1], [0, 0, 0, 0, 1, 1, 0, 1, 1], [0, 0, 0, 0, 1, 1, 1, 0, 1],
                    [0, 0, 0, 0, 1, 1, 1, 1, 1], [0, 0, 0, 1, 0, 0, 0, 1, 1], [0, 0, 0, 1, 0, 0, 1, 0, 1],
                    [0, 0, 0, 1, 0, 0, 1, 1, 1], [0, 0, 0, 1, 0, 1, 0, 0, 1], [0, 0, 0, 1, 0, 1, 0, 1, 1],
                    [0, 0, 0, 1, 0, 1, 1, 0, 1], [0, 0, 0, 1, 0, 1, 1, 1, 1], [0, 0, 0, 1, 1, 0, 0, 1, 1],
                    [0, 0, 0, 1, 1, 0, 1, 0, 1], [0, 0, 0, 1, 1, 0, 1, 1, 1], [0, 0, 0, 1, 1, 1, 0, 0, 1],
                    [0, 0, 0, 1, 1, 1, 0, 1, 1], [0, 0, 0, 1, 1, 1, 1, 0, 1], [0, 0, 0, 1, 1, 1, 1, 1, 1],
                    [0, 0, 1, 0, 0, 0, 0, 1, 1], [0, 0, 1, 0, 0, 0, 1, 0, 1], [0, 0, 1, 0, 0, 0, 1, 1, 1],
                    [0, 0, 1, 0, 0, 1, 0, 1, 1], [0, 0, 1, 0, 0, 1, 1, 0, 1], [0, 0, 1, 0, 0, 1, 1, 1, 1],
                    [0, 0, 1, 0, 1, 0, 0, 1, 1], [0, 0, 1, 0, 1, 0, 1, 0, 1], [0, 0, 1, 0, 1, 0, 1, 1, 1],
                    [0, 0, 1, 0, 1, 1, 0, 1, 1], [0, 0, 1, 0, 1, 1, 1, 0, 1], [0, 0, 1, 0, 1, 1, 1, 1, 1],
                    [0, 0, 1, 1, 0, 0, 1, 0, 1], [0, 0, 1, 1, 0, 0, 1, 1, 1], [0, 0, 1, 1, 0, 1, 0, 1, 1],
                    [0, 0, 1, 1, 0, 1, 1, 0, 1], [0, 0, 1, 1, 0, 1, 1, 1, 1], [0, 0, 1, 1, 1, 0, 1, 0, 1],
                    [0, 0, 1, 1, 1, 0, 1, 1, 1], [0, 0, 1, 1, 1, 1, 0, 1, 1], [0, 0, 1, 1, 1, 1, 1, 0, 1],
                    [0, 0, 1, 1, 1, 1, 1, 1, 1], [0, 1, 0, 0, 0, 0, 0, 1, 1], [0, 1, 0, 0, 0, 0, 1, 1, 1],
                    [0, 1, 0, 0, 0, 1, 0, 1, 1], [0, 1, 0, 0, 0, 1, 1, 1, 1], [0, 1, 0, 0, 1, 0, 0, 1, 1],
                    [0, 1, 0, 0, 1, 0, 1, 1, 1], [0, 1, 0, 0, 1, 1, 0, 1, 1], [0, 1, 0, 0, 1, 1, 1, 1, 1],
                    [0, 1, 0, 1, 0, 0, 0, 1, 1], [0, 1, 0, 1, 0, 0, 1, 1, 1], [0, 1, 0, 1, 0, 1, 0, 1, 1],
                    [0, 1, 0, 1, 0, 1, 1, 1, 1], [0, 1, 0, 1, 1, 0, 0, 1, 1], [0, 1, 0, 1, 1, 0, 1, 1, 1],
                    [0, 1, 0, 1, 1, 1, 0, 1, 1], [0, 1, 0, 1, 1, 1, 1, 1, 1], [0, 1, 1, 0, 0, 0, 1, 1, 1],
                    [0, 1, 1, 0, 0, 1, 1, 1, 1], [0, 1, 1, 0, 1, 0, 1, 1, 1], [0, 1, 1, 0, 1, 1, 1, 1, 1],
                    [0, 1, 1, 1, 0, 1, 1, 1, 1], [0, 1, 1, 1, 1, 1, 1, 1, 1], [1, 0, 0, 0, 0, 0, 0, 0, 0],
                    [1, 0, 0, 0, 1, 0, 0, 0, 0], [1, 0, 0, 1, 0, 0, 0, 0, 0], [1, 0, 0, 1, 0, 1, 0, 0, 0],
                    [1, 0, 0, 1, 1, 0, 0, 0, 0], [1, 0, 0, 1, 1, 1, 0, 0, 0], [1, 0, 1, 0, 0, 0, 0, 0, 0],
                    [1, 0, 1, 0, 0, 0, 1, 0, 0], [1, 0, 1, 0, 0, 1, 0, 0, 0], [1, 0, 1, 0, 0, 1, 1, 0, 0],
                    [1, 0, 1, 0, 1, 0, 0, 0, 0], [1, 0, 1, 0, 1, 0, 1, 0, 0], [1, 0, 1, 0, 1, 1, 0, 0, 0],
                    [1, 0, 1, 0, 1, 1, 1, 0, 0], [1, 0, 1, 1, 0, 0, 0, 0, 0], [1, 0, 1, 1, 0, 0, 1, 0, 0],
                    [1, 0, 1, 1, 0, 1, 0, 0, 0], [1, 0, 1, 1, 0, 1, 1, 0, 0], [1, 0, 1, 1, 1, 0, 0, 0, 0],
                    [1, 0, 1, 1, 1, 0, 1, 0, 0], [1, 0, 1, 1, 1, 1, 0, 0, 0], [1, 0, 1, 1, 1, 1, 1, 0, 0],
                    [1, 1, 0, 0, 0, 0, 0, 0, 0], [1, 1, 0, 0, 0, 0, 0, 1, 0], [1, 1, 0, 0, 0, 0, 1, 0, 0],
                    [1, 1, 0, 0, 0, 1, 0, 0, 0], [1, 1, 0, 0, 0, 1, 0, 1, 0], [1, 1, 0, 0, 1, 0, 0, 0, 0],
                    [1, 1, 0, 0, 1, 0, 0, 1, 0], [1, 1, 0, 0, 1, 0, 1, 0, 0], [1, 1, 0, 0, 1, 1, 0, 0, 0],
                    [1, 1, 0, 0, 1, 1, 0, 1, 0], [1, 1, 0, 1, 0, 0, 0, 0, 0], [1, 1, 0, 1, 0, 0, 0, 1, 0],
                    [1, 1, 0, 1, 0, 0, 1, 0, 0], [1, 1, 0, 1, 0, 1, 0, 0, 0], [1, 1, 0, 1, 0, 1, 0, 1, 0],
                    [1, 1, 0, 1, 0, 1, 1, 0, 0], [1, 1, 0, 1, 1, 0, 0, 0, 0], [1, 1, 0, 1, 1, 0, 0, 1, 0],
                    [1, 1, 0, 1, 1, 0, 1, 0, 0], [1, 1, 0, 1, 1, 1, 0, 0, 0], [1, 1, 0, 1, 1, 1, 0, 1, 0],
                    [1, 1, 0, 1, 1, 1, 1, 0, 0], [1, 1, 1, 0, 0, 0, 0, 0, 0], [1, 1, 1, 0, 0, 0, 0, 1, 0],
                    [1, 1, 1, 0, 0, 0, 1, 0, 0], [1, 1, 1, 0, 0, 0, 1, 1, 0], [1, 1, 1, 0, 0, 1, 0, 0, 0],
                    [1, 1, 1, 0, 0, 1, 0, 1, 0], [1, 1, 1, 0, 0, 1, 1, 0, 0], [1, 1, 1, 0, 1, 0, 0, 0, 0],
                    [1, 1, 1, 0, 1, 0, 0, 1, 0], [1, 1, 1, 0, 1, 0, 1, 0, 0], [1, 1, 1, 0, 1, 0, 1, 1, 0],
                    [1, 1, 1, 0, 1, 1, 0, 0, 0], [1, 1, 1, 0, 1, 1, 0, 1, 0], [1, 1, 1, 0, 1, 1, 1, 0, 0],
                    [1, 1, 1, 1, 0, 0, 0, 0, 0], [1, 1, 1, 1, 0, 0, 0, 1, 0], [1, 1, 1, 1, 0, 0, 1, 0, 0],
                    [1, 1, 1, 1, 0, 0, 1, 1, 0], [1, 1, 1, 1, 0, 1, 0, 0, 0], [1, 1, 1, 1, 0, 1, 0, 1, 0],
                    [1, 1, 1, 1, 0, 1, 1, 0, 0], [1, 1, 1, 1, 0, 1, 1, 1, 0], [1, 1, 1, 1, 1, 0, 0, 0, 0],
                    [1, 1, 1, 1, 1, 0, 0, 1, 0], [1, 1, 1, 1, 1, 0, 1, 0, 0], [1, 1, 1, 1, 1, 0, 1, 1, 0],
                    [1, 1, 1, 1, 1, 1, 0, 0, 0], [1, 1, 1, 1, 1, 1, 0, 1, 0], [1, 1, 1, 1, 1, 1, 1, 0, 0],
                    [1, 1, 1, 1, 1, 1, 1, 1, 0]]

    N = 8
    M = n // N

    mu = (M - m + 1) / (2 ** m)
    varWj = M * ((1.0 / (2 ** m)) - ((2.0 * m - 1.0) / (2.0 ** (2.0 * m))))

    sequence = choice(templates_m9)
    for jj in range(min(maxnumberOfTemplates, len(sequence))):
        j = 0
        Wj = [0] * N
        nu = [0] * (K + 1)
        for k in range(K + 1):
            nu[k] = 0
        for i in range(N):
            W_obs = 0
            for j in range(M - m + 1):
                match = 1
                for k in range(m):
                    if sequence[k] != int(bitstring[i * M + j + k]):
                        match = 0
                        break

                if match == 1:
                    W_obs += 1

            Wj[i] = W_obs

        j += m - 1

    chi2 = 0.0
    for i in range(N):
        chi2 += ((Wj[i] - mu) / (varWj ** 0.5)) ** 2
    p_value = scipy.special.gammaincc(N / 2.0, chi2 / 2.0)

    return p_value

def overlappingtemplatematching(bitstring):
    n = len(bitstring)

    #M = 20 for 1024 - 4096
    M = 35
    K = 5
    N = n//M
    m = 4

    pi = [0.0] * (6)

    print(f"MN = {M * N}\n n = {n}\n m = {m}\n N = {N}")

    sequence = [1] * m

    lam = (M-m+1)/(2 ** m)
    eta = lam / 2.0
    print(f"lam = {lam}\n 2 * lam = {2*lam}")

    def Pr(u, eta):
        if u == 0:
            p = math.exp(-eta)
        else:
            p = sum([math.exp(-eta - u * math.log(2) + x * math.log(eta) - math.log(scipy.special.gamma(x + 1))
                                  + math.log(scipy.special.gamma(u)) - math.log(scipy.special.gamma(x))
                                  - math.log(scipy.special.gamma(u - x + 1))) for x in range(1, u + 1)])
        return p

    pi_sum = 0.0
    for i in range(K):
        pi[i] = Pr(i, eta)
        pi_sum += pi[i]

    pi[K] = 1 - pi_sum

    print(f"pi = {pi} \n N * min(pi) = {N * min(pi)}")

    nu = [0] * (K+1)
    for i in range(N):
        W_obs = 0
        for j in range(int(M-m+1)):
            match = 1
            for k in range(m):
                if sequence[k] != int(bitstring[i * M + j + k]):
                    match = 0
            if match == 1:
                W_obs += 1
        if W_obs <= 4:
            nu[W_obs] += 1
        else:
            nu[K] += 1

    chi2 = 0.0
    for i in range(K+1):
        chi2 += ((nu[i] - N * pi[i]) ** 2) / (N * pi[i])

    p_value = scipy.special.gammaincc(K/2.0, chi2/2.0)

    return p_value

def universal(bitstring):
    n = len(bitstring)
    L = 2
    # Q = 4
    Q = 10 * 2 ** L
    K = n // L - Q


    expected_value = [0, 0.73264948, 1.5374383, 2.40160681, 3.31122472, 4.25342659, 5.2177052, 6.1962507, 7.1836656,
                 8.1764248, 9.1723243, 10.170032, 11.168765, 2.168070, 13.167693, 14.167488, 15.167379]
    variance = [0, 0.690, 1.338, 1.901, 2.358, 2.705, 2.954, 3.125, 3.238, 3.311, 3.356, 3.384, 3.401,
                3.410, 3.416, 3.419, 3.421]

    p = 2 ** L
    c = 0.7 - 0.8 / L + (4 + 32 / L) * K ** (-3 / L) / 15
    sigma = c * math.sqrt(variance[L]/K)
    sqrt2 = math.sqrt(2)
    phi_sum = 0.0
    T = [0] * (Q)
    for i in range(p):
        T[i] = 0
    for i in range(Q + 1):
        decRep = 0
        for j in range(L):
            decRep += int(bitstring[(i - 1) * L + j]) * 2 ** (L - 1 - j)
        T[decRep] = i


    i = Q + 1
    while i <= (Q + K):
        decRep = 0
        for j in range(L):
            decRep += int(bitstring[(i - 1) * L + j]) * 2 ** (L - 1 - j)
        phi_sum += math.log(i - T[decRep]) / math.log(2)
        T[decRep] = i
        i += 1

    phi = phi_sum / K
    arg = abs(phi - expected_value[L]) / (sqrt2 * sigma)
    p_value = math.erfc(arg)
    return p_value

def serial(bitstring):
    m = 2
    n = len(bitstring)

    def psi2(m, n):
        if m == 0 or m == -1:
            return 0.0
        nblocks = n
        powLen = 2 ** (m+1) - 1
        P = [0] * powLen

        i = 1
        while i < (powLen - 1):
            P[i] = 0
            i += 1


        for i in range(nblocks):
            k = 1
            for j in range(m):
                if int(bitstring[(i+j)%n]) == 0:
                    k *= 2
                elif int(bitstring[(i+j)%n]) == 1:
                    k = 2*k+1
            P[k-1] += 1

        psi2_sum = 0
        i = 2 ** m - 1
        while i < (2 ** (m+1) - 1):
            psi2_sum += P[i] ** 2
            i += 1
        psi2_sum = ((psi2_sum * 2 ** m) / n) - n
        return psi2_sum

    psim0 = psi2(m,n)
    psim1 = psi2(m-1,n)
    psim2 = psi2(m-2, n)
    del1 = psim0 - psim1
    del2 = psim0 - 2.0*psim1 + psim2
    p_value1 = scipy.special.gammaincc(2 ** (m-1)/2, del1/2.0)
    p_value2 = scipy.special.gammaincc(2 ** (m-2)/2, del2/2.0)
    return p_value1, p_value2

def approximate_entrophy(bitstring):
    #m = 1 for 64 and 128, 2 otherwise
    m = 2
    n = len(bitstring)

    ApEn = [0] * 2
    r = 0
    blocksize = m
    while blocksize <= (m+1):
        if blocksize == 0:
            ApEn[0] = 0.00
            r += 1
        else:
            nblocks = n
            powLen = 2 ** (blocksize + 1) - 1
            P = [0] * powLen
            i = 1
            while i < (powLen - 1):
                P[i] = 0
                i += 1
            for i in range(nblocks):
                k = 1
                for j in range(blocksize):
                    k <<= 1
                    if int(bitstring[(i+j)%n]) == 1:
                        k += 1
                P[k-1] += 1

            apen_sum = 0
            index = 2 ** blocksize - 1
            for i in range(2 ** blocksize):
                if P[index] > 0:
                    apen_sum += P[index] * math.log(P[index]/nblocks)
                index += 1

            apen_sum /= nblocks
            ApEn[r] = apen_sum
            r += 1
        blocksize += 1

    apen = ApEn[0] - ApEn[1]
    chi_2 = 2.0 * n * (math.log(2) - apen)
    p_value = scipy.special.gammaincc(2 ** (m-1), chi_2/2.0)
    return p_value

def cumulative_sums(bitstring):
    n = len(bitstring)
    S = 0
    sup = 0
    inf = 0

    def cephes_normal(x):
        if x > 0:
            arg = x / math.sqrt(2)
            result = 0.5 * (1 + math.erf(arg))
        else:
            arg = -x / math.sqrt(2)
            result = 0.5 * (1 - math.erf(arg))
        return result

    for k in range(n):
        S = S + 1 if int(bitstring[k]) else S - 1
        if S > sup:
            sup += 1
        if S < inf:
            inf -= 1
        z = sup if sup > -inf else -inf
        zrev = sup - S if sup - S > S - inf else S - inf


    sumf1 = 0.0
    k = (-n / z + 1) // 4
    while k <= ((n / z) - 1) // 4:
        sumf1 += cephes_normal(((4 * k + 1)*z)/math.sqrt(n))
        sumf1 -= cephes_normal(((4 * k - 1)*z)/math.sqrt(n))
        k += 1

    sumf2 = 0.0
    k  = (-n / z - 3) // 4
    while k <= ((n/z)-1) // 4:
        sumf2 += cephes_normal(((4 * k+ 3)* z)/math.sqrt(n))
        sumf2 -= cephes_normal(((4 * k + 1)* z)/math.sqrt(n))
        k += 1

    sumb1 = 0.0
    k = (-n / zrev + 1) // 4
    while k <= (n / zrev - 1) // 4:
        sumb1 += cephes_normal(((4 * k + 1) * zrev) / math.sqrt(n))
        sumb1 -= cephes_normal(((4 * k - 1) * zrev) / math.sqrt(n))
        k += 1

    sumb2 = 0.0
    k = (-n / zrev - 3) // 4
    while k <= (n / zrev - 1) // 4:
        sumb2 += cephes_normal(((4 * k + 3) * zrev) / math.sqrt(n))
        sumb2 -= cephes_normal(((4 * k + 1) * zrev) / math.sqrt(n))
        k += 1

    p_valuef = 1.0 - sumf1 + sumf2
    p_valueb = 1.0 - sumb1 + sumb2
    return p_valuef, p_valueb

def random_excursions(bitstring):
    n = len(bitstring)
    J = 0
    S_k = [0] * n
    cycle = [0] * n
    nu = numpy.ndarray(shape=(6, 8), dtype=int)
    pi = [[0.0000000000, 0.00000000000, 0.00000000000, 0.00000000000, 0.00000000000, 0.0000000000],
          [0.5000000000, 0.25000000000, 0.12500000000, 0.06250000000, 0.03125000000, 0.0312500000],
          [0.7500000000, 0.06250000000, 0.04687500000, 0.03515625000, 0.02636718750, 0.0791015625],
          [0.8333333333, 0.02777777778, 0.02314814815, 0.01929012346, 0.01607510288, 0.0803755143],
          [0.8750000000, 0.01562500000, 0.01367187500, 0.01196289063, 0.01046752930, 0.0732727051]]
    stateX = [-4, -3, -2, -1, 1, 2, 3, 4]
    counter = [0, 0, 0, 0, 0, 0, 0, 0]

    S_k[0] = 2 * int(bitstring[0]) - 1
    i = 1
    while i < n:
        S_k[i] = S_k[i - 1] + 2 * int(bitstring[i]) - 1
        if S_k[i] == 0:
            J += 1
            cycle[J] = i
        i += 1

    if S_k[n - 1] != 0:
        J += 1
    cycle[J] = n

    cycleStart = 0
    cycleStop = cycle[1]
    for k in range(6):
        for i in range(8):
            nu[k, i] = 0
    j = 1
    while j <= J:
        for i in range(8):
            counter[i] = 0
        i = cycleStart
        while i < cycleStop:
            if (S_k[i] >= 1 and S_k[i] <= 4) or (S_k[i] >= -4 and S_k[i] <= -1):
                if S_k[i] < 0:
                    b = 4
                else:
                    b = 3
                counter[S_k[i] + b] += 1
            i += 1

        cycleStart = cycle[j] + 1
        if j < J:
            cycleStop = cycle[j + 1]

        for i in range(8):
            if counter[i] >= 0 and counter[i] <= 4:
                nu[counter[i], i] += 1
            elif counter[i] >= 5:
                nu[5, i] += 1
        j += 1

    p_list = [0.0] * 8
    for i in range(8):
        x = stateX[i]
        chi2 = 0.0
        for k in range(6):
            chi2 += (nu[k, i] - J * pi[abs(x)][k]) ** 2 / (J * pi[abs(x)][k])
        p_list[i] = scipy.special.gammaincc(2.5, chi2/2.0)

    p_value_sum = 0.0
    for x in p_list:
        if x < 0.01:
            p_value_sum = 0.0
            break
        else:
            p_value_sum += x

    return p_value_sum

def random_excursions_variant(bitstring):
    n = len(bitstring)
    stateX = [-9, -8, -7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    J = 0
    S_k = [0] * n
    S_k[0] = 2 * int(bitstring[0]) - 1
    i = 1
    while i < n:
        S_k[i] = S_k[i - 1] + 2 * int(bitstring[i]) - 1
        if S_k[i] == 0:
            J += 1
        i += 1

    if S_k[n - 1] != 0:
        J += 1

    p_list = [0] * 18
    for p in range(18):
        x = stateX[p]
        count = 0
        for i in range(n):
            if S_k[i] == x:
                count += 1
        p_list[p] = math.erfc(abs(count-J)/(math.sqrt(2.0 * J * (4.0 * abs(x) - 2))))

    p_value_sum = 0.0
    for x in p_list:
        if x < 0.01:
            p_value_sum = 0.0
            break
        else:
            p_value_sum += x


    return p_value_sum

if __name__ == "__main__":
    t = "0110110101"
    e = "1100100100001111110110101010001000100001011010001100001000110100110001001100011001100010100010111000"
    e_128 = "1100110000010101011011000100110011100000000000100100110101010001000100111101011010" \
            "0000001101011111001100111001101101100010110010"
    e_50 = "10111011110010110100011100101110111110000101101001"

    with open("data.e", "r") as e1:
        e_binary_rank = "".join([line.strip() for line in e1.readlines()])

    # Tested and work for 64 and above:
    # print(f"Frequency test: {freq_monobit(e)}")
    # print(f"Block test: {block_freq(10, e)}")
    # print(f"Runs Test: {runs(e)}")
    # print(f"Longest Run Test: {longestrunofones(e_128)}")
    # print(f"Discrete Fourier Transform: {discrete_fourier_transform(e_binary_rank[:256])}")
    # print(f"NonOverlappingTemplate Matching: {nonOverlappingTemplateMatchings(e_binary_rank[:256])}")
    # print(f'Serial: {serial(e_binary_rank[:64])}')
    # print(f'Cusum: {cumulative_sums(e_binary_rank[:64])}')
    # print(f'RE: {random_excursions(e_binary_rank[:256])}')
    # print(f'REV: {random_excursions_variant(e_binary_rank[:64])}')

    # Do not work for under 128:
    # print(f'Approximate Entrophy: {approximate_entrophy(e_binary_rank[:128])}')

    # Do not work for under 256:
    # print(f"Universal: {universal(e_binary_rank[:256])}")
    # print(f"Rank test: {binary_matrix_rank(e_binary_rank[:256])}")

    # Do not work for under 1024:
    # print(f"OverlappingTemplate Matching: {overlappingtemplatematching(e_binary_rank[:4096])}")
