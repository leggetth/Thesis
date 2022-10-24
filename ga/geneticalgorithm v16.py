from geneticalgorithm import geneticalgorithm as ga
from nistsp80022tests import *
import pickle
import timer


def f(ga_array):
    freq_monobit_p_value = freq_monobit(ga_array)
    block_freq_p_value = block_freq(ga_array)
    runs_p_value = runs(ga_array)
    long_run_p_value = longestrunofones(ga_array)
    bin_mat_ran_p_value = binary_matrix_rank(ga_array)
    dis_four_tran_p_value = discrete_fourier_transform(ga_array)
    nonover_temp_mat_p_value = nonOverlappingTemplateMatchings(ga_array)
    uni_p_value = universal(ga_array)
    ser_p_value1, ser_p_value2 = serial(ga_array)
    app_en_p_value = approximate_entrophy(ga_array)
    cu_sum_p_value1, cu_sum_p_value2 = cumulative_sums(ga_array)
    re_p_value = random_excursions(ga_array)
    rev_p_value = random_excursions_variant(ga_array)
    tot_p_value = freq_monobit_p_value + block_freq_p_value + runs_p_value + long_run_p_value + bin_mat_ran_p_value + \
                  dis_four_tran_p_value + nonover_temp_mat_p_value + uni_p_value + ser_p_value1 + ser_p_value2 + \
                  app_en_p_value + cu_sum_p_value1 + cu_sum_p_value2 + re_p_value + rev_p_value
    return -tot_p_value

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
    dimension=2048,
    variable_type='bool',
    algorithm_parameters=algorithm_param,
    convergence_curve=True
)

t = timer.Timer()
t.start()
model.run()
print(t.stop())

# output_dict = {}
# time_lst = []
# t = timer.Timer()
# for x in range(1, 11):
#     t.start()
#     model.run()
#     time_lst.append(t.stop())
#     output = "".join(str(x) for x in model.output_dict["variable"].astype("int"))
#     output_dict[x] = output
#
# with open("freq_monobit_bin.txt", "wb") as f1:
#     pickle.dump([output_dict, time_lst], f1)
#
# #print("\n\n", [x for x in output.astype("int")])