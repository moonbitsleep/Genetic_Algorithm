"""
多种群遗传算法
添加粗略局部搜索
"""
import copy
import math
import random
from operator import attrgetter

import numpy as np

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
POP_SIZE = 100  # 种群大小
TOTAL_POP = 10  # 总种群个数10个，9个子种群，1个最优保存种群
GROUP = 30  # 选择权重，百分之百
CROSS_RATE = 0.8
MUTATE_RATE = 0.2
GEN_MAX = 20  # 连续GEN_MAX代保持之前最优值


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
            assert len(records) == project_num  # 确保每个人都完成了8个项目
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
            (parent1, parent2) = self.members[i], self.members[i + 1]
            child1, child2 = self.crossover(parent1, parent2)
            self.mutate(child1)
            self.mutate(child2)
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

    def get_weak(self):
        return max(self.members, key=attrgetter('fitness'))

    def get_best(self):
        self.get_fitness()
        return min(self.members, key=attrgetter('fitness'))

    def select(self):
        if self.members[-1].fitness is None:
            self.get_fitness()  # 选择前计算适应度
        num_to_select = math.floor(self.size * (GROUP / 100))
        sample = random.sample(range(self.size), num_to_select)
        sample_members = sorted([self.members[i] for i in sample], key=attrgetter('fitness'))
        return sample_members[:2]

    def crossover(self, parent1, parent2):
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
    def mutate(chromosome):
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

    def get_fitness(self):
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
        """粗略局部搜索生成新dna"""
        # pidx = random.randint(1, total_people)
        p_projects = [(pidx - 1) * project_num + i for i in range(1, project_num + 1)]
        # 拿到索引和项目
        p_idxs, p_seqs = self.get_proj_idx_seq(seq_copy, p_projects)
        random.shuffle(p_seqs)
        self.change_seq(seq_copy, p_idxs, p_seqs)
        return Chromosome(seq_copy)

def MGA():

    best_fitness = 99999999
    best_dna_seq = None
    metrics = None
    gen0 = 0
    gen = 0
    pops = [Population(POP_SIZE) for _ in range(TOTAL_POP)]

    while gen0 <= GEN_MAX:
        gen = gen + 1  # 记录迭代次数

        # 子种群内部进化
        for i in range(TOTAL_POP - 1):
            pop = pops[i]
            pop.evolve_population()

        # 迁移进化
        for j in range(TOTAL_POP - 1):
            pop1 = pops[j]
            pop2 = pops[(j + 1) % (TOTAL_POP - 1)]
            parent1 = pop1.get_best()
            for i in range(POP_SIZE):
                parent2 = pop2.members[i]
                child1, child2 = pop2.crossover(parent1, parent2)
                pop2.mutate(child1)
                pop2.mutate(child2)
                parent1.compute_fitness()
                parent2.compute_fitness()
                child1.compute_fitness()
                child2.compute_fitness()
                better = min([parent2, child1, child2], key=attrgetter('fitness'))
                idx = pop2.members.index(parent2)
                pop2.members[idx] = better

        # 优秀dna填充到最优保存种群
        elite_pop = pops[-1]
        elite_pop.get_fitness()

        best_dnas = []
        for i in range(TOTAL_POP - 1):
            best_dnas.append(pops[i].get_best())
        for dna in best_dnas:
            weak = elite_pop.get_weak()
            if weak.fitness > dna.fitness:
                elite_pop.members.remove(weak)
                elite_pop.members.append(dna)
                # print("elite_pop size:", len(elite_pop.members))

        # 计算最终最优解
        # 局部搜索
        for mem in elite_pop.members:
            elite_pop.local_search(mem)
        final_best_dna = elite_pop.get_best()
        print("Gen: {}, best dna seq: {}".format(gen, final_best_dna.sequence))
        print("best fitness: {}".format(final_best_dna.fitness))
        print("makespan: {}, avg_total_wait: {}, T_W_wait: {}".format(final_best_dna.makespan,
                                                                      final_best_dna.total_wait / total_people,
                                                                      final_best_dna.greater_than_threshold))
        # 记录最优解
        if final_best_dna.fitness < best_fitness:
            best_fitness = final_best_dna.fitness
            best_dna_seq = final_best_dna.sequence
            metrics = [final_best_dna.makespan, final_best_dna.total_wait / total_people,
                       final_best_dna.greater_than_threshold]
            gen0 = 0  # 保持次数清0
            print("最优解:", best_fitness)
            print("最优序列: ", best_dna_seq)
            print("最优指标: ", metrics)
        else:
            gen0 = gen0 + 1
            print("已经重复最优值{}次".format(gen0))


    print("===========================================")
    print("最优解:", best_fitness)
    print("最优序列: ", best_dna_seq)
    print("最优指标: ", metrics)
    print("===========================================")

if __name__ == '__main__':
    MGA()