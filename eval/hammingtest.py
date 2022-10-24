import pickle
import numpy
import csv

def hamming_dist(x, y):
    assert len(x) == len(y)

    hamming = 0
    for i in range(len(x)):
        char_xor = ord(x[i]) ^ ord(y[i])
        hamming += bin(char_xor).count("1")

    return hamming

def hamming_dist_byte(x, y):
    assert len(x) == len(y)

    hamming = 0
    for i in range(len(x)):
        char_xor = int(x[i]) ^ int(y[i])
        hamming += bin(char_xor).count("1")

    return hamming

def hamming_test_binary(input):
    ham_list_1 = []
    ham_list_2 = []
    ham_list_3 = []
    ham_list_4 = []
    ham_list_5 = []
    ham_list_6 = []
    ham_list_7 = []
    ham_list_8 = []
    ham_list_9 = []
    ham_list_10 = []

    for x in range(1, 11):
        if x == 1:
            for i in range(1, 11):
                ham_list_1.append(hamming_dist(input.get(x), input.get(i)))

        if x == 2:
            for i in range(1, 11):
                ham_list_2.append(hamming_dist(input.get(x), input.get(i)))

        if x == 3:
            for i in range(1, 11):
                ham_list_3.append(hamming_dist(input.get(x), input.get(i)))

        if x == 4:
            for i in range(1, 11):
                ham_list_4.append(hamming_dist(input.get(x), input.get(i)))

        if x == 5:
            for i in range(1, 11):
                ham_list_5.append(hamming_dist(input.get(x), input.get(i)))

        if x == 6:
            for i in range(1, 11):
                ham_list_6.append(hamming_dist(input.get(x), input.get(i)))

        if x == 7:
            for i in range(1, 11):
                ham_list_7.append(hamming_dist(input.get(x), input.get(i)))

        if x == 8:
            for i in range(1, 11):
                ham_list_8.append(hamming_dist(input.get(x), input.get(i)))

        if x == 9:
            for i in range(1, 11):
                ham_list_9.append(hamming_dist(input.get(x), input.get(i)))

        if x == 10:
            for i in range(1, 11):
                ham_list_10.append(hamming_dist(input.get(x), input.get(i)))

    sum_list = [ham_list_1, ham_list_2, ham_list_3, ham_list_4, ham_list_5, ham_list_6,
                ham_list_7, ham_list_8, ham_list_9, ham_list_10]

    avg_list = []
    for x in sum_list:
        average = sum(i for i in x)
        avg_list.append(average)

    return sum(avg_list) / 100


    # return {"list": sum_list,
    #         "average": sum(sum_list) / len(sum_list),
    #         "avg_byte": sum([(x / 8) for x in sum_list]) / len(sum_list)}

def hamming_test_bytes(input):
    sum_list = []
    for x in range(1, 10):
        if x == 1:
            for i in range(2, 11):
                sum_list.append(hamming_dist_byte(input.get(x), input.get(i)))

        if x == 2:
            for i in range(3, 11):
                sum_list.append(hamming_dist_byte(input.get(x), input.get(i)))

        if x == 3:
            for i in range(4, 11):
                sum_list.append(hamming_dist_byte(input.get(x), input.get(i)))

        if x == 4:
            for i in range(5, 11):
                sum_list.append(hamming_dist_byte(input.get(x), input.get(i)))

        if x == 5:
            for i in range(6, 11):
                sum_list.append(hamming_dist_byte(input.get(x), input.get(i)))

        if x == 6:
            for i in range(7, 11):
                sum_list.append(hamming_dist_byte(input.get(x), input.get(i)))

        if x == 7:
            for i in range(8, 11):
                sum_list.append(hamming_dist_byte(input.get(x), input.get(i)))

        if x == 8:
            for i in range(9, 11):
                sum_list.append(hamming_dist_byte(input.get(x), input.get(i)))

        if x == 9:
            for i in range(10, 11):
                sum_list.append(hamming_dist_byte(input.get(x), input.get(i)))

    return {"list": sum_list,
            "average": sum(sum_list) / len(sum_list),
            "avg_byte": sum([(x / 8) for x in sum_list]) / len(sum_list)}

def avg_pvalue(x):
    return sum(x.values()) / len(x)

if __name__ == "__main__":
    filename = "test.csv"


    bin_dict = {1: "101100", 2: "110011", 3: "101010", 4: "001100", 5: "110011",
                6: "001010", 7: "110101", 8: "010101", 9: "100111", 10: "111000"}

    # with open("freq_monobit_bin.txt", "rb") as f1:
    #     output_dict_freq, output_dict_p, time_lst_freq = pickle.load(f1)
    # # print(output_dict,"\n", time_lst)
    #
    # # ham_test_freq = hamming_test_binary(output_dict_freq)
    # # print(f"Frequency List: {ham_test_freq['list']}")
    # # print(f"Frequency Average: {ham_test_freq['average']}")
    # # print(f"Frequency Average Byte: {ham_test_freq['avg_byte']}")
    #
    # # print(hamming_test_binary(output_dict_freq))
    #
    # # with open("test.csv", "w") as csvfile:
    # #     csvwriter = csv.writer(csvfile)
    # #     csvwriter.writerow(hamming_test_binary(output_dict_freq))
    # #     csvwriter.writerow(time_lst_freq)
    #
    # print(hamming_test_binary(output_dict_freq))
    # print(avg_pvalue(output_dict_p))
    # print(f"Frequency Time Average: {sum([float(x) for x in time_lst_freq]) / 10}")
    #
    # with open("freq_block_bin.txt", "rb") as f2:
    #     output_dict_fb, output_dict_pb, time_lst_fb = pickle.load(f2)
    # # print(output_dict_2, "\n",time_lst_2)
    #
    # # ham_test_freq_block = hamming_test_binary(output_dict_fb)
    # # print(f"\nFrequency Block List: {ham_test_freq_block['list']}")
    # # print(f"Frequency Block Average: {ham_test_freq_block['average']}")
    # # print(f"Frequency Block Average Byte: {ham_test_freq_block['avg_byte']}")
    #
    # print(hamming_test_binary(output_dict_fb))
    # print(avg_pvalue(output_dict_pb))
    # print(f"Frequency Block Time Average: {sum([float(x) for x in time_lst_fb]) / 10}")

    # with open("test.csv", "w") as csvfile:
    #     csvwriter = csv.writer(csvfile)
    #     csvwriter.writerow(hamming_test_binary(output_dict_fb))
    #     csvwriter.writerow(time_lst_fb)

    # with open("runs_bin.txt", "rb") as f3:
    #     output_dict_runs, output_dict_rp, time_lst_runs = pickle.load(f3)
    #
    # print(hamming_test_binary(output_dict_runs))
    # print(avg_pvalue(output_dict_rp))
    # print(f"Runs Time Average: {sum([float(x) for x in time_lst_runs]) / 10}")
    #
    # with open("lroo_bin.txt", "rb") as f4:
    #     output_dict_lroo, output_dict_lroop, time_lst_lroo = pickle.load(f4)
    #
    # print(hamming_test_binary(output_dict_lroo))
    # print(avg_pvalue(output_dict_lroop))
    # print(f"LROO Time Average: {sum([float(x) for x in time_lst_lroo]) / 10}")
    #
    # with open("bmr_bin.txt", "rb") as f5:
    #     output_dict_bmr, output_dict_bmrp, time_lst_bmr = pickle.load(f5)
    #
    # print(hamming_test_binary(output_dict_bmr))
    # print(avg_pvalue(output_dict_bmrp))
    # print(f"BMR Time Average: {sum([float(x) for x in time_lst_bmr]) / 10}")
    #
    # with open("dft_bin.txt", "rb") as f6:
    #     output_dict_dft, output_dict_dftp, time_lst_dft = pickle.load(f6)
    #
    # print(hamming_test_binary(output_dict_dft))
    # print(avg_pvalue(output_dict_dftp))
    # print(f"DFT Time Average: {sum([float(x) for x in time_lst_dft]) / 10}")
    #
    # with open("notm_bin.txt", "rb") as f7:
    #     output_dict_notm, output_dict_notmp, time_lst_notm = pickle.load(f7)
    #
    # print(hamming_test_binary(output_dict_notm))
    # print(avg_pvalue(output_dict_notmp))
    # print(f"NOTM Time Average: {sum([float(x) for x in time_lst_notm]) / 10}")
    #
    # with open("universal_bin.txt", "rb") as f8:
    #     output_dict_universal, output_dict_universalp, time_lst_universal = pickle.load(f8)
    #
    # print(hamming_test_binary(output_dict_universal))
    # print(avg_pvalue(output_dict_universalp))
    # print(f"Universal Time Average: {sum([float(x) for x in time_lst_universal]) / 10}")

    with open("serial_bin.txt", "rb") as f9:
        output_dict_serial, output_dict_serialp, time_lst_serial = pickle.load(f9)

    print(hamming_test_binary(output_dict_serial))
    print(avg_pvalue(output_dict_serialp))
    print(f"Serial Time Average: {sum([float(x) for x in time_lst_serial]) / 10}")

    with open("aent_bin.txt", "rb") as f10:
        output_dict_aent, output_dict_aentp, time_lst_aent = pickle.load(f10)

    print(hamming_test_binary(output_dict_aent))
    print(avg_pvalue(output_dict_aentp))
    print(f"AENT Time Average: {sum([float(x) for x in time_lst_aent]) / 10}")

    with open("cusum_bin.txt", "rb") as f11:
        output_dict_cusum, output_dict_cusump, time_lst_cusum = pickle.load(f11)

    print(hamming_test_binary(output_dict_cusum))
    print(avg_pvalue(output_dict_cusump))
    print(f"Cusum Time Average: {sum([float(x) for x in time_lst_cusum]) / 10}")

    with open("ranexc_bin.txt", "rb") as f12:
        output_dict_ranexc, output_dict_ranexcp, time_lst_ranexc = pickle.load(f12)

    print(hamming_test_binary(output_dict_ranexc))
    print(avg_pvalue(output_dict_ranexcp))
    print(f"Ranexc Time Average: {sum([float(x) for x in time_lst_ranexc]) / 10}")

    with open("ranexcvar_bin.txt", "rb") as f13:
        output_dict_ranexcvar, output_dict_ranexcvarp, time_lst_ranexcvar = pickle.load(f13)

    print(hamming_test_binary(output_dict_ranexcvar))
    print(avg_pvalue(output_dict_ranexcvarp))
    print(f"Ranexcvar Time Average: {sum([float(x) for x in time_lst_ranexcvar]) / 10}")


