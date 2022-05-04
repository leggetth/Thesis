from geneticalgorithm import geneticalgorithm as ga
import math
import nistmatrixdefinitions as nmatrix
import pickle
import timer


def binary_matrix_rank(ga_array):
    n = len(ga_array)
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

        p_30 = 1 - (p_32 + p_31)

        F_32 = 0
        F_31 = 0
        for k in range(N):
            matrix = nmatrix.BinaryMatrix(ga_array, M, Q, k)
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
    return -p_value


algorithm_param = {
    'max_num_iteration': 500,
    'population_size': 50,
    'mutation_probability': 0.05,
    'elit_ratio': 0.02,
    'crossover_probability': 0.25,
    'parents_portion': 0.15,
    'crossover_type': 'two_point',
    'max_iteration_without_improv': None}


model=ga(
    function=binary_matrix_rank,
    dimension=128,
    variable_type='bool',
    algorithm_parameters=algorithm_param,
    convergence_curve=False
)

# model.run()
output_dict = {}
time_lst = []
t = timer.Timer()
for x in range(1, 11):
    t.start()
    model.run()
    time_lst.append(t.stop())
    output = "".join(str(x) for x in model.output_dict["variable"].astype("int"))
    output_dict[x] = output

with open("bmr_bin.txt", "wb") as f1:
    pickle.dump([output_dict, time_lst], f1)
#