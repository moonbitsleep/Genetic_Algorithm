import copy
import time
from operator import attrgetter

from common_setting import *
from pop import Population
from collections import deque
import matplotlib.pyplot as plt

class MultiPopGeneticAlgorithm:
    def __init__(self, pop_num):
        self.pop_num = pop_num
        self.save_pop = Population(POP_SIZE, CROSS_RATE, MUTATE_RATE)
        self.best_fit = 99999
        self.best_dna = None
        self.generation_fits = deque()
        self.pops = deque()
        self.init_pops()

    def init_pops(self):
        for i in range(self.pop_num - 1):
            pop = Population(POP_SIZE, CROSS_RATE, MUTATE_RATE)
            self.pops.append(pop)

    def run(self):
        gen0 = 0
        gen = 0

        while gen0 <= GEN_MAX:
            gen = gen + 1

            # self evolve
            for pop in self.pops:
                pop.evolve()

            # migration evolution
            for j in range(len(self.pops)):
                pop1 = self.pops[j]
                pop2 = self.pops[(j + 1) % (self.pop_num - 1)]
                parent1 = pop1.get_best_chrom()
                for i in range(POP_SIZE):
                    parent2 = pop2.members[i]
                    child1, child2 = pop2.crossover(parent1, parent2)
                    pop2.mutate(child1)
                    pop2.mutate(child2)
                    parent1.compute_metric()
                    parent2.compute_metric()
                    child1.compute_metric()
                    child2.compute_metric()
                    better = min([parent2, child1, child2], key=attrgetter('fitness'))
                    idx = pop2.members.index(parent2)
                    pop2.members[idx] = better

            #save
            self.save_pop.get_fitness()
            best_dnas = deque()
            for i in range(len(self.pops)):
                best_dnas.append(self.pops[i].get_best_chrom())

            for dna in best_dnas:
                weak = self.save_pop.get_weak_chrom()
                if weak.fitness > dna.fitness:
                    self.save_pop.members.remove(weak)
                    self.save_pop.members.append(dna)

            # final best
            final_best_dna = self.save_pop.get_best_chrom()
            print("===========================")
            print("Gen:", gen)
            final_best_dna.print_metric()
            print("===========================")
            self.generation_fits.append(final_best_dna.fitness)
            if self.best_fit > final_best_dna.fitness:
                self.best_fit = final_best_dna.fitness
                self.best_dna = final_best_dna
                gen0 = 0
            else:
                gen0 = gen0 + 1
                print("repeat {} times".format(gen0))

    def print_final_best(self):
        print("*********best*************")
        self.best_dna.print_metric()
        print("*********best*************")

    def draw_plot(self):
        plt.plot(self.generation_fits)
        plt.show()


if __name__ == '__main__':
    mga = MultiPopGeneticAlgorithm(10)
    start = time.time()
    mga.run()
    end = time.time()
    print("waste time:", end - start)
    mga.print_final_best()
    mga.draw_plot()
