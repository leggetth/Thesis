from nistsp80022tests import runs, cumulative_sums, approximate_entrophy, universal
from geneticalgorithm import geneticalgorithm as ga
from hammingtest import hamming_dist
from numpy import arange
import timer

class Sensitivity_Analysis():
    def __init__(self, fitness_functions, w1, w2, w3, w4):
        self.fitness_functions = fitness_functions
        self.w1 = w1
        self.w2 = w2
        self.w3 = w3
        self.w4 = w4

    def f(self, ga_array):
        p_value_list = []
        for fitness_function in self.fitness_functions:
            p_value_list.append(fitness_function(ga_array))

        total_p_value = (self.w1 * p_value_list[0]) + (self.w2 * (p_value_list[1][0] + p_value_list[1][1])) + \
                        (self.w3 * p_value_list[2]) + (self.w4 * p_value_list[3])

        return -total_p_value

    def genetic_algorithm(self):
        algorithm_param = {
            'max_num_iteration': 500,
            'population_size': 17,
            'mutation_probability': 0.05,
            'elit_ratio': 1 / 17,
            'crossover_probability': 0.25,
            'parents_portion': 3 / 17,
            'crossover_type': 'two_point',
            'max_iteration_without_improv': None}

        model = ga(
            function=self.f,
            dimension=256,
            variable_type='bool',
            algorithm_parameters=algorithm_param,
            convergence_curve=False
        )

        t = timer.Timer()

        t.start()
        model.run()
        output_time = t.stop()

        output_key = "".join(str(x) for x in model.output_dict["variable"].astype("int"))

        output_p_value = model.output_dict["function"]

        return {"key": output_key, "time": output_time, "p_value": output_p_value}

    def evaluation(self):
        keys = []
        times = []
        p_values = []
        for i in range(2):
            ga_output = self.genetic_algorithm()
            keys.append(ga_output["key"])
            times.append(ga_output["time"])
            p_values.append(ga_output["p_value"])

        hamming_distance = hamming_dist(keys[0], keys[1])
        avg_time = sum([float(x) for x in times]) / 2
        avg_p_value = sum(p_values) / 2

        eval_value = hamming_distance + (1/avg_time)

        return {"eval_value": eval_value, "hamming_distance": hamming_distance,
                "avg_time": avg_time, "avg_p_value": avg_p_value}

if __name__ == "__main__":

    best_eval_value = 0.0
    for w1 in arange(0.0, 4/3, 1/3):
        for w2 in arange(0.0, 4/3, 1/3):
            for w3 in arange(0.0, 4/3, 1/3):
                for w4 in arange(0.0, 4/3, 1/3):
                    sa = Sensitivity_Analysis([runs, cumulative_sums, approximate_entrophy, universal],
                                              w1, w2, w3, w4)
                    eval_output = sa.evaluation()
                    if eval_output["eval_value"] > best_eval_value:
                        best_weights = {"w1": w1, "w2": w2, "w3": w3, "w4": w4}
                        best_hamming_distance = eval_output["hamming_distance"]
                        best_avg_time = eval_output["avg_time"]
                        best_avg_p_value = eval_output["avg_p_value"]
                        best_eval_value = eval_output["eval_value"]

    print(f"\n\nBest Values:\nWeights: {best_weights}\nHamming Distance: {best_hamming_distance}\n"
          f"Time: {best_avg_time}\nP_value: {best_avg_p_value}")