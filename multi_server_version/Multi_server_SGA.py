"""
假设每个科室有2个服务台
每人依然检查8个项目
顾客同时到达
服务时长固定
染色体编码不变，解码不采用半活动解码，优先选择人少的队列插入记录
相当于只考虑调度顺序的变化带来的后果，不考虑具体队列的情况
"""
import copy
import math
import random
from operator import attrgetter

import numpy as np
from matplotlib import pyplot as plt

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
resource_look_up = [
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
]
total_people = 60
project_num = len(cost_time_lookup)
T_W = 15  # 等待阈值
GAMMA1 = 1  # 惩罚系数1，用于实验1
GAMMA2 = 10  # 惩罚系数2，用于超阈值等待
POP_SIZE = 150  # 种群大小
GROUP = 30  # 选择权重，百分之百
CROSS_RATE = 0.8
MUTATE_RATE = 1.0
N_GENERATIONS = 500


class Chromosome:
    def __init__(self, sequence):
        self.sequence = sequence  # list DNA序列
        self.fitness = None  # 适应度值
        self.makespan = None  # 完工时间
        self.total_wait = None  # 总等待时间
        self.greater_than_threshold = None  # 超阈值等待时间
        self.service_idle = None  # 服务台空闲时间
        self.patient_table = [[] for _ in range(total_people)]  # 顾客存储[project_index, service_index, start, end]
        self.service_table = [
            [[] for _ in range(resource_look_up[sidx])] for sidx in range(project_num)
        ]  # 服务台存储[patient_index, start, end]

    @staticmethod
    def translate_operation(opt):
        """解码操作"""
        opt = opt - 1
        return opt // project_num, opt % project_num

    @staticmethod
    def get_people_last_end(people_records):
        """获取顾客什么时候有空的时间"""
        if len(people_records) == 0:
            return 0
        else:
            people_records.sort(key=lambda e: e[3])
            return people_records[-1][3]

    @staticmethod
    def get_early_queue(queues_records):
        """
        找到服务台最早空闲队列
        :param queues_records: 服务台队列
        :return: [队列索引, 最早空闲时间]
        """
        targets = [0, 99999]
        for i in range(len(queues_records)):
            queue = queues_records[i]
            if len(queue) == 0:
                return i, 0
            else:
                queue.sort(key=lambda e: e[2])
                if queue[-1][2] < targets[1]:
                    targets[0] = i
                    targets[1] = queue[-1][2]
        return targets[0], targets[1]

    def translate(self):
        """
        解码
        :return: None
        """
        for operation in self.sequence:
            patient_index, project_index = self.translate_operation(operation)
            cost_time = cost_time_lookup[project_index]
            patient_records = self.patient_table[patient_index]  # 客户的所有记录
            queues_records = self.service_table[project_index]  # 项目服务台多个队列的所有记录

            # 获取客户最早有空时间
            patient_early_time = self.get_people_last_end(patient_records)
            # 获取最早空闲的队列及时间
            queue_index, project_early_time = self.get_early_queue(queues_records)
            # 确定占用时间片段
            start_time = max(patient_early_time, project_early_time)
            end_time = start_time + cost_time
            # 添加记录
            patient_records.append([project_index, queue_index, start_time, end_time])
            queues_records[queue_index].append([patient_index, start_time, end_time])

    def clean_table_cache(self):
        self.patient_table = [[] for _ in range(total_people)]  # 顾客存储[project_index, service_index, start, end]
        self.service_table = [
            [[] for _ in range(resource_look_up[sidx])] for sidx in range(project_num)
        ]  # 服务台存储[patient_index, start, end]

    def compute_fitness(self):
        """
        解码
        :return: None
        """
        if self.fitness is not None:
            self.clean_table_cache()
        self.translate()
        W_sum = 0
        w_thanT_sum = 0
        tmp_lates = []
        for records in self.patient_table:
            assert len(records) == project_num # 确保每个人都完成了8个项目
            records.sort(key=lambda e: e[3])
            tmp_lates.append(records[-1][3])
            for i in range(len(records)):
                if i == 0:
                    wait_time = records[i][2] - 0
                    if wait_time == 0:
                        continue
                    W_sum += wait_time
                    if wait_time - T_W > 0:
                        w_thanT_sum += (wait_time - T_W)
                else:
                    wait_time = records[i][2] - records[i - 1][3]
                    W_sum += wait_time
                    if wait_time - T_W > 0:
                        w_thanT_sum += (wait_time - T_W)
        maxF = max(tmp_lates)
        self.makespan = maxF
        self.total_wait = W_sum
        self.greater_than_threshold = w_thanT_sum
        self.fitness = maxF + W_sum + GAMMA2 * w_thanT_sum

    def print_metric(self):
        """打印指标"""
        print("seq: ", self.sequence)
        print("obj: ", self.fitness)
        print("makespan: ", self.makespan)
        print("total_wait: ", self.total_wait)
        print("greater_than_threshold: ", self.greater_than_threshold)

    def print_status(self):
        """打印表格"""
        print("service records:")
        for service in self.service_table:
            for records in service:
                print(records)
        print("patient records:")
        for records in self.patient_table:
            print(records)

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


    def kill_weak(self):
        weakest = max(self.members, key=attrgetter('fitness'))
        self.members.remove(weakest)

    def get_best(self):
        self._get_fitness()
        return sorted(self.members, key=attrgetter('fitness'))[:1]


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
        global project_num
        global MUTATE_RATE
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


if __name__ == '__main__':
    best_fitness = 99999999
    population = Population(POP_SIZE)
    best_fits = []
    metrics = None
    for i in range(N_GENERATIONS):
        population.evolve_population()
        dna = population.get_best()
        best = dna[0]
        print("Gen: {}, best dna seq: {}".format(i, best.sequence))
        print("best fitness: {}".format(best.fitness))
        print("makespan: {}, total_wait: {}, T_W_wait: {}".format(best.makespan, best.total_wait,
                                                                  best.greater_than_threshold))
        if best_fitness > best.fitness:
            best_fitness = best.fitness
            metrics = [best.makespan, best.total_wait, best.greater_than_threshold]
        best_fits.append(best_fitness)

    plt.plot(best_fits)
    plt.show()
    print(best_fitness)
    print(metrics)

