import copy
import numpy

class BinaryMatrix():
    def __init__(self, bitstring, M, Q, k):
        self.rows = M
        self.columns = Q
        self.matrix = self.def_matrix(bitstring, k)
        self.base_rank = min(M, Q)

    def perform_row_operations(self, i, forward_elimination):
        if forward_elimination == 1:
            j = i + 1
            while j < self.rows:
                if self.matrix[j, i] == 1:
                    for k in range(self.columns):
                        self.matrix[j, k] = (self.matrix[j, k] + self.matrix[i, k]) % 2
                j += 1

        else:
            j = i - 1
            while j >= 0:
                if self.matrix[j, i] == 1:
                    for k in range(self.columns):
                        self.matrix[j, k] = (self.matrix[j, k] + self.matrix[i, k]) % 2
                j -= 1

    def find_unit_element_and_swap(self, i, forward_elimination):
        row_op = 0
        if forward_elimination == 1:
            index = i + 1
            while index < self.rows and self.matrix[index, i] == 0:
                index += 1
            if index < self.rows:
                row_op = self.swap_rows(i, index)
        else:

            index = i - 1
            while index >= 0 and self.matrix[index, i] == 0:
                index -= 1
            if index >= 0:
                row_op = self.swap_rows(i, index)

        return row_op

    def swap_rows(self, i, index):
        temp = copy.copy(self.matrix[i, :])
        self.matrix[i, :] = self.matrix[index, :]
        self.matrix[index, :] = temp
        return 1

    def determine_rank(self):
        rank = self.base_rank
        for i in range(self.rows):
            allzeroes = 1
            for j in range(self.columns):
                if self.matrix[i, j] == 1:
                    allzeroes = 0
                    break

            if allzeroes == 1:
                rank -= 1
        return rank

    def compute_rank(self):
        i = 0
        while i < (self.base_rank - 1):
            if self.matrix[i, i] == 1:
                self.perform_row_operations(i, 1)
            else:
                if self.find_unit_element_and_swap(i, 1) == 1:
                    self.perform_row_operations(i, 1)
            i += 1

        i = self.base_rank - 1
        while i > 0:
            if self.matrix[i, i] == 1:
                self.perform_row_operations(i, 0)
            else:
                if self.find_unit_element_and_swap(i, 0) == 1:
                    self.perform_row_operations(i, 0)
            i -= 1

        return self.determine_rank()

    def def_matrix(self, bitstring, k):
        matrix = numpy.ndarray(shape=(self.rows, self.columns), dtype=int)
        for i in range(self.rows):
            for j in range(self.columns):
                matrix[i, j] = int(bitstring[k * (self.rows * self.columns) + j + i * self.rows])
        return matrix