from geneticalgorithm import geneticalgorithm as ga
import scipy.special
import pickle
import timer


def block_freq(ga_array, block_len=32):
    ga_array = "".join(str(x) for x in ga_array.astype("int"))
    num_block = len(ga_array) // block_len
    sum = 0
    for i in range(num_block):
        block_sum = 0
        for j in range(num_block):
            block_sum += int(ga_array[j+i*num_block])
        pi = block_sum / num_block
        v = pi - 0.5
        sum += v * v
    chi_squared = 4.0 * block_len * sum
    p_value = scipy.special.gammaincc((num_block/2.0), (chi_squared/2.0))
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
    function=block_freq,
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

with open("freq_block_bin.txt", "wb") as f1:
    pickle.dump([output_dict, time_lst], f1)