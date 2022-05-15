import copy
import random
import time

import numpy as np

from common_setting import *
from collections import deque
from chrom import Chromosome
from operator import attrgetter

class Population:
    def __init__(self, pop_size, cross_rate, mutate_rate):
        self.pop_size = pop_size
        self.pc = cross_rate
        self.pm = mutate_rate
        self.members = deque(maxlen=pop_size)
        self._init_members1()

    def _init_members1(self):
        """random init"""
        for i in range(self.pop_size):
            sequence = deque()
            normal_bound = [i + 1 for i in range(BATCH_PEOPLE * PROJECT_NUM)]
            end_people = TOTAL_PEOPLE - (BATCH_COUNT - 1) * BATCH_PEOPLE
            end_bound = [i + 1 for i in range(end_people * PROJECT_NUM)]
            for i in range(BATCH_COUNT - 1):
                batch = deque(np.random.choice(normal_bound, size=BATCH_PEOPLE * PROJECT_NUM, replace=False))
                # print(len(batch))
                sequence.append(batch)
            end_seq = deque(np.random.choice(end_bound, size=end_people * PROJECT_NUM, replace=False))
            # print(len(end_seq))
            # for i in range((BATCH_PEOPLE - end_people) * PROJECT_NUM):
            #     end_seq.append(0)
            # print(len(end_seq))
            sequence.append(end_seq)
            # print(sequence)
            self.members.append(Chromosome(sequence))

    def _init_members2(self):
        """crossover and mutate"""
        pass

    def get_fitness(self):
        for mem in self.members:
            mem.compute_metric()

    @staticmethod
    def get_proj_idx_seq(sequence, targets):
        idxs = deque()
        orders = deque()
        for i in range(len(sequence)):
            if sequence[i] in targets:
                idxs.append(i)
                orders.append(sequence[i])
        return idxs, orders

    @staticmethod
    def change_seq(child, change_idxs, new_seqs):
        for i in range(len(child)):
            if i in change_idxs:
                child[i] = new_seqs.popleft()

    def crossover(self, parent1, parent2):
        if random.random() < self.pc:
            child1_seq = copy.deepcopy(parent1.seq)
            child2_seq = copy.deepcopy(parent2.seq)

            batch_idx = random.randint(0, BATCH_COUNT - 1)
            # TODO: add param for chromosome last batch size
            pidx = random.randint(1, BATCH_PEOPLE - 1)
            p_projects = [(pidx - 1) * PROJECT_NUM + i for i in range(1, PROJECT_NUM + 1)]

            p1_idxs, p1_seqs = self.get_proj_idx_seq(parent1.seq, p_projects)
            p2_idxs, p2_seqs = self.get_proj_idx_seq(parent2.seq, p_projects)
            self.change_seq(child1_seq[batch_idx], p1_idxs, p2_seqs)
            self.change_seq(child2_seq[batch_idx], p2_idxs, p1_seqs)
            return Chromosome(child1_seq), Chromosome(child2_seq)
        else:
            return parent1, parent2

    def mutate(self, chromosome):
        if random.random() < self.pm:
            batch_idx = random.randint(0, BATCH_COUNT-1)
            batch_seq = chromosome.seq[batch_idx]
            idx1, idx2 = np.random.choice(len(batch_seq), size=2, replace=False)
            batch_seq[idx1], batch_seq[idx2] = batch_seq[idx2], batch_seq[idx1]

    def get_best_chrom(self):
        self.get_fitness()
        return min(self.members, key=attrgetter("fitness"))

    def evolve(self):
        # TODO: select
        for i in range(self.pop_size - 1):
            (parent1, parent2) = self.members[i], self.members[i + 1]
            child1, child2 = self.crossover(parent1, parent2)
            self.mutate(child1)
            self.mutate(child2)
            parent1.compute_metric()
            parent2.compute_metric()
            child1.compute_metric()
            child2.compute_metric()
            better = min([parent1, child1, child2], key=attrgetter('fitness'))
            idx = self.members.index(parent1)
            self.members[idx] = better

def test_crossover():
    s = 0
    Pop = Population(POP_SIZE, CROSS_RATE, MUTATE_RATE)
    for i in range(10000):
        idx = random.randint(0, Pop.pop_size - 2)
        parent1 = Pop.members[idx]
        parent2 = Pop.members[idx + 1]
        child1, child2 = Pop.crossover(parent1, parent2)
        parent1.compute_metric()
        parent2.compute_metric()
        child1.compute_metric()
        child2.compute_metric()
        bef = max([parent1, parent2], key=attrgetter("fitness"))
        aft = min([child1, child2], key=attrgetter("fitness"))
        print(bef.fitness, aft.fitness)
        if bef.fitness > aft.fitness:
            s += 1
    print("调优成功{}次".format(s))

def test_evolve():
    Pop = Population(POP_SIZE, CROSS_RATE, MUTATE_RATE)
    print(Pop.get_best_chrom().fitness)
    Pop.evolve()
    print(Pop.get_best_chrom().fitness)

if __name__ == '__main__':
    s1 = time.time()
    test_evolve()
    e1 = time.time()
    print("迭代一轮花费{}s".format(e1 - s1))

