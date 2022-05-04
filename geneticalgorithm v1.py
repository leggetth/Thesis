from geneticalgorithm import geneticalgorithm as ga
import math
import pickle
import timer


def freq_monobit(ga_array):
    array = [int(x) for x in ga_array]
    ones = array.count(1)
    zeros = array.count(0)
    difference = abs(ones - zeros)
    s_obs = float(difference) / math.sqrt(float(len(ga_array)))
    p_value = math.erfc(s_obs / math.sqrt(2))
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
    function=freq_monobit,
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

with open("freq_monobit_bin.txt", "wb") as f1:
    pickle.dump([output_dict, time_lst], f1)
#
#print("\n\n", [x for x in output.astype("int")])