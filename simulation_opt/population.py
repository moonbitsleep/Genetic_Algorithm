import copy
import math
from operator import attrgetter

import numpy as np

from params import *
import random
from chrom import Chromosome

class Population:
    def __init__(self, size):
        self.size = size
        self.members = []
        self.init_population()

    def init_population(self):
        """初始化种群"""
        sequence_size = total_people * project_num
        for i in range(self.size):
            sequence = random.sample(range(1, sequence_size + 1), sequence_size)
            self.members.append(Chromosome(sequence))
            # 局部搜索
        for mem in self.members:
            self.local_search(mem)

    def evolve_population(self):
        """种群迭代"""
        global POP_SIZE
        for i in range(POP_SIZE - 1):
            (parent1, parent2) = self.members[i], self.members[i+1]
            child1, child2 = self._crossover(parent1, parent2)
            self._mutate(child1)
            self._mutate(child2)
            parent1.compute_fitness()
            parent2.compute_fitness()
            child1.compute_fitness()
            child2.compute_fitness()
            better = min([parent1, child1, child2], key=attrgetter('fitness'))
            idx = self.members.index(parent1)
            self.members[idx] = better
        # 局部搜索
        for mem in self.members:
            self.local_search(mem)

    def kill_weak(self):
        weakest = max(self.members, key=attrgetter('fitness'))
        self.members.remove(weakest)

    def get_best(self):
        self._get_fitness()
        return min(self.members, key=attrgetter('fitness'))


    def select(self):
        if self.members[-1].fitness is None:
            self._get_fitness()  # 选择前计算适应度
        num_to_select = math.floor(self.size * (GROUP/100))
        sample = random.sample(range(self.size), num_to_select)
        sample_members = sorted([self.members[i] for i in sample], key=attrgetter('fitness'))
        return sample_members[:2]

    def _crossover(self, parent1, parent2):
        """得到两个child_seq, 再构造染色体返回"""
        # 选择一个顾客id
        if random.random() < CROSS_RATE:
            pidx = random.randint(1, total_people)
            p_projects = [(pidx - 1) * project_num + i for i in range(1, project_num + 1)]
            child1_seq = copy.deepcopy(parent1.sequence)
            child2_seq = copy.deepcopy(parent2.sequence)
            # 拿到parent1上该顾客的所有项目及顺序
            p1_idxs, p1_seqs = self.get_proj_idx_seq(parent1.sequence, p_projects)
            # 拿到parent2上该顾客的所有项目及顺序
            p2_idxs, p2_seqs = self.get_proj_idx_seq(parent2.sequence, p_projects)
            random.shuffle(p1_seqs)
            random.shuffle(p2_seqs)
            self.change_seq(child1_seq, p1_idxs, p2_seqs)
            self.change_seq(child2_seq, p2_idxs, p1_seqs)
            return Chromosome(child1_seq), Chromosome(child2_seq)
        else:
            return parent1, parent2

    @staticmethod
    def _mutate(chromosome):
        """变异"""
        if random.random() < MUTATE_RATE:
            idx1, idx2 = np.random.choice(len(chromosome.sequence), size=2, replace=False)
            chromosome.sequence[idx1], chromosome.sequence[idx2] = chromosome.sequence[idx2], chromosome.sequence[idx1]

    @staticmethod
    def change_seq(child, change_idxs, new_seqs):
        for i in range(len(child)):
            if i in change_idxs:
                child[i] = new_seqs.pop(0)

    @staticmethod
    def get_proj_idx_seq(sequence, targets):
        idxs = []
        orders = []
        for i in range(len(sequence)):
            if sequence[i] in targets:
                idxs.append(i)
                orders.append(sequence[i])
        return idxs, orders

    def _get_fitness(self):
        """统一计算整个种群的fitness"""
        for mem in self.members:
            mem.compute_fitness()

    @staticmethod
    def get_large_waits_pidx(chromosome):
        """从染色体中得到较大等待顾客的pid"""
        if chromosome.fitness is None:
            chromosome.compute_fitness()
        people_records = chromosome.patient_table
        pidxs = []
        for pid in range(len(people_records)):
            records = people_records[pid]
            records.sort(key=lambda e: e[3])
            for i in range(len(records)):
                if i == 0:
                    wait_time = records[i][2] - 0
                    if wait_time == 0:
                        continue
                    if wait_time >= 20:
                        pidxs.append(pid)
                        break
                else:
                    wait_time = records[i][2] - records[i - 1][3]
                    if wait_time - T_W >= 20:
                        pidxs.append(pid)
                        break
        return pidxs

    # def local_search(self, chromosome):
    #     """局部搜索"""
    #     pidxs = self.get_large_waits_pidx(chromosome)
    #     tmp1 = self.generate_new_dna(copy.deepcopy(chromosome.sequence), random.choice(pidxs))
    #     tmp2 = self.generate_new_dna(copy.deepcopy(chromosome.sequence), random.choice(pidxs))
    #     tmp3 = self.generate_new_dna(copy.deepcopy(chromosome.sequence), random.choice(pidxs))
    #     tmp4 = self.generate_new_dna(copy.deepcopy(chromosome.sequence), random.choice(pidxs))
    #     tmp1.compute_fitness()
    #     tmp2.compute_fitness()
    #     tmp3.compute_fitness()
    #     tmp4.compute_fitness()
    #     chromosome.compute_fitness()
    #     better = min([tmp1, tmp2, tmp3, tmp4, chromosome], key=attrgetter('fitness'))
    #     idx = self.members.index(chromosome)
    #     self.members[idx] = better

    def local_search(self, chromosome):
        neighbor1 = chromosome.explore_search1()
        neighbor1.compute_fitness()
        chromosome.compute_fitness()
        better = min([neighbor1, chromosome], key=attrgetter('fitness'))
        idx = self.members.index(chromosome)
        self.members[idx] = better

    def generate_new_dna(self, seq_copy, pidx):
        # pidx = random.randint(1, total_people)
        p_projects = [(pidx - 1) * project_num + i for i in range(1, project_num + 1)]
        # 拿到索引和项目
        p_idxs, p_seqs = self.get_proj_idx_seq(seq_copy, p_projects)
        random.shuffle(p_seqs)
        self.change_seq(seq_copy, p_idxs, p_seqs)
        return Chromosome(seq_copy)