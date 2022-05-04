from geneticalgorithm import geneticalgorithm as ga
import math
import pickle
import timer


def runs(ga_array):
    n = len(ga_array)
    array = [int(x) for x in ga_array]
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
    function=runs,
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

with open("runs_bin.txt", "wb") as f1:
    pickle.dump([output_dict, time_lst], f1)
#