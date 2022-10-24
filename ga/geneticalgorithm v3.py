from geneticalgorithm import geneticalgorithm as ga
import pickle
import timer
from nistsp80022tests import runs

def f(ga_array):
    runs_p_value = runs(ga_array)
    return -runs_p_value


algorithm_param = {
    'max_num_iteration': 500,
    'population_size': 17,
    'mutation_probability': 0.05,
    'elit_ratio': 1/17,
    'crossover_probability': 0.25,
    'parents_portion': 3/17,
    'crossover_type': 'two_point',
    'max_iteration_without_improv': None}


model=ga(
    function=f,
    dimension=256,
    variable_type='bool',
    algorithm_parameters=algorithm_param,
    convergence_curve=False
)

# model.run()
output_dict_key = {}
output_dict_p_value= {}
time_lst = []
t = timer.Timer()
for x in range(1, 11):
    t.start()
    model.run()
    time_lst.append(t.stop())
    output_key = "".join(str(x) for x in model.output_dict["variable"].astype("int"))
    output_dict_key[x] = output_key
    output_dict_p_value[x] = model.output_dict["function"]

with open("runs_bin.txt", "wb") as f1:
    pickle.dump([output_dict_key, output_dict_p_value, time_lst], f1)
#