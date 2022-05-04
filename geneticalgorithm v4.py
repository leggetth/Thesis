from geneticalgorithm import geneticalgorithm as ga
import scipy.special
import pickle
import timer


def longestrunofones(ga_array):
    n = len(ga_array)
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
            if int(ga_array[i*M+j]) == 1:
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
    function=longestrunofones,
    dimension=64,
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

with open("lroo_bin.txt", "wb") as f1:
    pickle.dump([output_dict, time_lst], f1)
#