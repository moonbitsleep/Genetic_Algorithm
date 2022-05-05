"""
1. 随机产生初始解，建立初始种群，种群规模为Z；进行解码和适应度计算
2. 使用交叉算子和突变算子对当前种群中对解进行重组，产生后代，形成新的种群。进行解码和适应度计算
3. 迭代：重复1-2，直到解一定代数内最优值没有发生明显改进，输出当前最优解
"""
import math
import time
import matplotlib.pyplot as plt
import random
from operator import attrgetter
import copy
import numpy as np
from A1B1 import generate_dna

cost_time_lookup = [
    3,  # 0体质测试   4
    3,  # 1内科      4
    4,  # 2外科      3
    2,  # 3眼耳口鼻科 5
    3,  # 4验血      4
    2,  # 5心电图    5
    5,  # 6X光      3
    6,  # 7B超      2
]  # 共28分钟

# cost_time_lookup = [
#     3, 3, 4, 6
# ]

total_people = 40
project_num = len(cost_time_lookup)
T_W = 15  # 等待阈值
GAMMA = 5 # 惩罚系数
POP_SIZE = 150  # 种群大小
GROUP = 30  # 选择权重，百分之百
CROSS_RATE = 0.8
MUTATE_RATE = 1.0
N_GENERATIONS = 500

def translate_operation(opt):
    """解码操作"""
    opt = opt - 1
    return opt // project_num, opt % project_num


class Chromosome:
    def __init__(self, sequence):
        self.sequence = sequence  # list DNA序列
        self.fitness = None
        self.makespan = None
        self.total_wait = None
        self.greater_than_threshold = None

    def translate(self):
        global cost_time_lookup
        global total_people
        global project_num
        people_records = [[] for _ in range(total_people)]
        project_records = [[] for _ in range(project_num)]
        for operation in self.sequence:
            people_index, project_index = translate_operation(operation)
            cost_time = cost_time_lookup[project_index]

            people = people_records[people_index]  # 客户的所有记录
            project = project_records[project_index]  # 项目的所有记录
            # 创建记录
            people_last_end_time = self.get_people_last_end(people)  # 人有空的最早时间
            project_last_ends = self.get_project_last_end(project)

            # 科室第一条记录
            if type(project_last_ends) == int:
                start_time = people_last_end_time
                end_time = start_time + cost_time
                people.append([project_index, start_time, end_time])
                project.append([people_index, start_time, end_time])

            if type(project_last_ends) == list:

                tmp_start = self.get_middle_start_time(project_last_ends, cost_time, people_last_end_time)
                if tmp_start == -1:
                    # 科室末尾插入
                    start_time = max(people_last_end_time, project_last_ends[-1])
                    end_time = start_time + cost_time
                    people.append([project_index, start_time, end_time])
                    project.append([people_index, start_time, end_time])
                else:
                    # 科室中间插入
                    start_time = tmp_start
                    end_time = start_time + cost_time
                    people.append([project_index, start_time, end_time])
                    project.append([people_index, start_time, end_time])
        return people_records, project_records

    def compute_fitness(self):
        global GAMMA
        global T_W
        people_records, project_records = self.translate()
        W_sum = 0
        w_thanT_sum = 0
        tmp_lates = []
        for records in people_records:
            records.sort(key=lambda e:e[2])
            tmp_lates.append(records[-1][2])
            for i in range(len(records)):
                if i == 0:
                    wait_time = records[i][1] - 0
                    if wait_time == 0:
                        continue
                    W_sum += wait_time
                    if wait_time - T_W > 0:
                        w_thanT_sum += (wait_time - T_W)
                else:
                    wait_time = records[i][1] - records[i - 1][2]
                    W_sum += wait_time
                    if wait_time - T_W > 0:
                        w_thanT_sum += (wait_time - T_W)
        maxF = max(tmp_lates)
        self.makespan = maxF
        self.total_wait = W_sum
        self.greater_than_threshold = w_thanT_sum
        self.fitness = maxF + W_sum + GAMMA * w_thanT_sum
        # print("****************")
        # print("dna seq:", self.sequence)
        # print("people's records: ")
        # for p in people_records:
        #     print(p)
        # print("service records: ")
        # for s in project_records:
        #     s.sort(key=lambda e: e[2])
        #     print(s)
        # print("fitness:", self.fitness)
        # print("makespan:", maxF)
        # print("total_wait:", W_sum)
        # print("greater_than_threshold:", w_thanT_sum)
        # print("****************")

    @staticmethod
    def get_people_last_end(people_table):
        if len(people_table) == 0:
            return 0
        else:
            people_table.sort(key=lambda e: e[2])
            return people_table[-1][2]

    @staticmethod
    def get_project_last_end(project_table):
        """得到项目表中每个人的结束时间"""
        if len(project_table) == 0:
            return 0
        else:
            project_table.sort(key=lambda e: e[2])
            return [record[-1] for record in project_table]

    @staticmethod
    def get_middle_start_time(end_time_list, cost_time, people_start):
        """
        如果二者间时间间隔大于2倍的检查时间，表明有足够的时间
        如果后者 - 客户开始时间仍大于2倍的检查时间(要减去之前分配的一次检查)，则表明可以插入
        否则，只能在末尾插入
        """
        for i in range(1, len(end_time_list)):
            if end_time_list[i] - end_time_list[i - 1] >= 2 * cost_time:
                if end_time_list[i] - people_start >= 2 * cost_time:
                    return max(people_start, end_time_list[i - 1])
        return -1


class Population:
    def __init__(self, size):
        self.size = size
        self.members = []
        self.__seed_population()

    # def __seed_population(self):
    #     global total_people
    #     global project_num
    #     sequence_size = total_people * project_num
    #     sample_seq = [i for i in range(1, sequence_size + 1)]
    #     init_seq = self.get_init_seq()
    #
    #     init_len = len(init_seq)
    #     for i in init_seq:
    #         sample_seq.remove(i)
    #     random.shuffle(init_seq)
    #     for i in range(self.size):
    #         sequence = copy.deepcopy(init_seq)
    #         rest = random.sample(sample_seq, sequence_size - init_len)
    #         sequence.extend(rest)
    #         self.members.append(Chromosome(sequence))

    def __seed_population(self):
        global total_people
        global project_num
        sequence_size = total_people * project_num
        for i in range(self.size):
            sequence = random.sample(range(1, sequence_size + 1), sequence_size)
            self.members.append(Chromosome(sequence))
            # sequence = generate_dna()
            # self.members.append(Chromosome(sequence))
        # 局部搜索
        for mem in self.members:
            self.local_search(mem)

    # def get_init_seq(self):
    #     actions = [[j + i * project_num for i in range(total_people)] for j in range(1, project_num + 1)]
    #     res = []
    #     res.extend(actions[0][:4])
    #     res.extend(actions[1][4:8])
    #     res.extend(actions[2][8:11])
    #     res.extend(actions[3][11:16])
    #     res.extend(actions[4][16:20])
    #     res.extend(actions[5][20:25])
    #     res.extend(actions[6][25:28])
    #     res.extend(actions[7][28:])
    #     assert len(res) == 30
    #     return res

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
        global total_people
        global project_num
        global CROSS_RATE
        if random.random() < CROSS_RATE:
            pidx = random.randint(1, total_people)
            p_projects = [(pidx - 1) * project_num + i for i in range(1, project_num + 1)]
            child1_seq = copy.deepcopy(parent1.sequence)
            child2_seq = copy.deepcopy(parent2.sequence)
            # 拿到parent1上该顾客的所有项目及顺序
            p1_idxs, p1_seqs = self.get_proj_idx_seq(parent1.sequence, p_projects)
            # 拿到parent2上该顾客的所有项目及顺序
            p2_idxs, p2_seqs = self.get_proj_idx_seq(parent2.sequence, p_projects)

            # random.shuffle(p1_seqs)
            # random.shuffle(p2_seqs)

            self.change_seq(child1_seq, p1_idxs, p2_seqs)
            self.change_seq(child2_seq, p2_idxs, p1_seqs)
            return Chromosome(child1_seq), Chromosome(child2_seq)
        else:
            return parent1, parent2

    @staticmethod
    def _mutate(chromosome):
        """变异"""
        global project_num
        global MUTATE_RATE
        if random.random() < MUTATE_RATE:
            idx1, idx2 = np.random.choice(len(chromosome.sequence), size=2, replace=False)
            chromosome.sequence[idx1], chromosome.sequence[idx2] = chromosome.sequence[idx2], chromosome.sequence[idx1]


    def local_search(self, chromosome):
        """局部搜索"""
        pidxs = self.get_large_waits_pidx(chromosome)
        tmp1 = self.generate_new_dna(copy.deepcopy(chromosome.sequence), random.choice(pidxs))
        tmp2 = self.generate_new_dna(copy.deepcopy(chromosome.sequence), random.choice(pidxs))
        tmp3 = self.generate_new_dna(copy.deepcopy(chromosome.sequence), random.choice(pidxs))
        tmp4 = self.generate_new_dna(copy.deepcopy(chromosome.sequence), random.choice(pidxs))
        tmp1.compute_fitness()
        tmp2.compute_fitness()
        tmp3.compute_fitness()
        tmp4.compute_fitness()
        chromosome.compute_fitness()
        better = min([tmp1, tmp2, tmp3, tmp4, chromosome], key=attrgetter('fitness'))
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

    @staticmethod
    def get_large_waits_pidx(chromosome):
        """从染色体中得到较大等待顾客的pid"""
        people_records, project_records = chromosome.translate()
        pidxs = []
        for pid in range(len(people_records)):
            records = people_records[pid]
            records.sort(key=lambda e:e[2])
            for i in range(len(records)):
                if i == 0:
                    wait_time = records[i][1] - 0
                    if wait_time == 0:
                        continue
                    if wait_time >= 20:
                        pidxs.append(pid)
                        break
                else:
                    wait_time = records[i][1] - records[i - 1][2]
                    if wait_time - T_W >= 20:
                        pidxs.append(pid)
                        break
        return pidxs


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

if __name__ == '__main__':
    best_fitness = 99999999
    population = Population(POP_SIZE)
    best_fits = []
    metrics = None
    for i in range(N_GENERATIONS):
        population.evolve_population()
        best = population.get_best()
        print("Gen: {}, best dna seq: {}".format(i, best.sequence))
        print("best fitness: {}".format(best.fitness))
        print("makespan: {}, total_wait: {}, T_W_wait: {}".format(best.makespan, best.total_wait, best.greater_than_threshold))
        if best_fitness > best.fitness:
            best_fitness = best.fitness
            metrics = [best.makespan, best.total_wait, best.greater_than_threshold]
        best_fits.append(best_fitness)

    plt.plot(best_fits)
    plt.show()
    print(best_fitness)
    print(metrics)