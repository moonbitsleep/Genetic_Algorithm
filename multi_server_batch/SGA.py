import time

from common_setting import *
from pop import Population
from collections import deque

class SinglePopGeneticAlgorithm:
    def __init__(self):
        self.pop = Population(POP_SIZE, CROSS_RATE, MUTATE_RATE)
        self.best_fit = 99999
        self.generation_fits = deque()
        self.best_chrom = None

    def run(self, steps):
        for i in range(steps):
            self.pop.evolve()
            best_dna = self.pop.get_best_chrom()
            print("===========================")
            print("Gen:", i)
            best_dna.print_metric()
            print("===========================")
            self.generation_fits.append(best_dna.fitness)
            if best_dna.fitness < self.best_fit:
                self.best_fit = best_dna.fitness
                self.best_chrom = best_dna

    def print_final_best(self):
        print("*********best*************")
        self.best_chrom.print_metric()
        print("*********best*************")

if __name__ == '__main__':
    algo = SinglePopGeneticAlgorithm()
    s1 = time.time()
    algo.run(150)
    e1 = time.time()
    print("waste time:", e1 - s1)
    algo.print_final_best()
