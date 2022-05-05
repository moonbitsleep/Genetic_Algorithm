"""
单种群遗传算法
染色体chromosome:
"""
import time

import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt


TIME_MAX = 1000

class PatientRecords:
    """客户必须+1，因为0不是合法下标"""

    def __init__(self, total_num):
        self.total = total_num
        self.store = self.__init_store()

    def __init_store(self):
        # 舍弃0位置,不保存
        store = [[] for _ in range(self.total + 1)]
        return store

    def add_record(self, pidx, record):
        self.check_self(pidx, record)
        self.store[pidx].append(record)

    def check_self(self, pidx, record):
        st = record[1]
        records = self.store[pidx]
        for record in records:
            if record[1] == st:
                print("p-insert-record:", record)
                print("p-records:", records)
                raise ValueError("同一时间在两个科室接受服务")

    def get_last_project_end(self, pidx):
        """得到上一项目的结束时间"""
        records = self.store[pidx]
        if len(records) == 0:
            return 0
        records.sort(key=lambda e: e[2])
        return records[-1][2]

    def not_full(self):
        for records in self.store:
            if len(records) == 0:
                return True
        return False


class ServiceRecords:
    """服务台数目+1，因为0号不是合法服务台"""

    def __init__(self, total_num, serve_times):
        """
        :param total_num: 科室数目
        :param serve_times: 一个列表，保存平均服务时间; 或者一个矩阵，保存所有客户所有项目的服务时间
        """
        self.total = total_num
        self.serve_times = serve_times
        self.store = self.__init_store()

    def __init_store(self):
        # 舍弃0位置,不保存
        store = [[] for _ in range(self.total)]
        return store

    def not_full(self):
        for records in self.store:
            if len(records) == 0:
                return True
        return False

    def add_record(self, service_idx, record):
        """
        :param service_idx: 服务台编号
        :param record: [客户index, 开始检查时间, 结束检查时间]
        """
        self.check_self(service_idx, record)
        self.store[service_idx].append(record)

    def check_self(self, service_idx, record):
        st = record[1]
        records = self.store[service_idx]
        for record in records:
            if st == record[1]:
                print("s-insert-record:", record)
                print("s-records:", records)
                raise ValueError("同一时刻服务一个客户")

    def get_idle_times(self, service_idx):
        """得到该服务台的空闲时间区间"""
        service = self.store[service_idx]
        idles = []
        if len(service) == 0:
            idles.append([0, TIME_MAX])
            return idles
        service.sort(key=lambda e: e[1])  # 按照开始检查时间排序
        # 检查0时刻和第一个病人到达时是否有空闲——一般没有
        if service[0][1] - 0 > 0:
            idles.append([0, service[0][1]])
        # 中间的空闲区间
        for i in range(0, len(service)):
            if i == len(service) - 1:
                break
            second_time = service[i + 1][1]
            first_time = service[i][2]
            if second_time - first_time > 0:
                idles.append([first_time, second_time])
        # 末尾的空闲时间
        # idles.append([service[-1][2], TIME_MAX])

        return idles


class DNA:
    def __init__(self, people_num, service_times):
        self.people_num = people_num
        self.service_times = service_times
        self.bound = people_num * (len(service_times) - 1)
        self.patient_records = PatientRecords(people_num)
        self.service_records = ServiceRecords(len(service_times), service_times)
        self.fitness = 0  # 适应度值
        self.genes = self.__init_dna()

    def __init_dna(self):
        return np.random.choice([i for i in range(1, self.bound + 1)], size=self.bound, replace=False)

    def __getitem__(self, item):
        """返回dna保存的值"""
        return self.genes[item]

    def __repr__(self):
        genes = list(self.genes)
        return genes.__str__()

    def clear_records(self):
        self.patient_records = PatientRecords(self.people_num)
        self.service_records = ServiceRecords(len(self.service_times), self.service_times)


def compute_fitness(dna):
    """
    计算一条DNA解的适应度值
    顾客总等待时间 + max(科室结束服务时间) + 超阈值等待时间
    """
    patient_status = dna.patient_records.store
    service_status = dna.service_records.store
    max_F = get_max_makespan(service_status)
    T_W = 5
    gamma = 5
    W_sum = 0
    w_thanT_sum = 0
    for records in patient_status:
        for i in range(0, len(records)):
            if i == 0:
                wait_time = records[i][1] - 0
                W_sum += wait_time
                if wait_time - T_W > 0:
                    w_thanT_sum += wait_time - T_W
            else:
                wait_time = records[i][1] - records[i - 1][2]
                W_sum += wait_time
                if wait_time - T_W > 0:
                    w_thanT_sum += wait_time - T_W
    fitn = max_F + W_sum + w_thanT_sum * gamma
    return fitn


def compute_metric(dna):
    patient_status = dna.patient_records.store
    service_status = dna.service_records.store
    max_F = get_max_makespan(service_status)
    T_W = 15
    gamma = 5
    W_sum = 0
    w_thanT_sum = 0
    for records in patient_status:
        for i in range(0, len(records)):
            if i == 0:
                wait_time = records[i][1] - 0
                W_sum += wait_time
                if wait_time - T_W > 0:
                    w_thanT_sum += wait_time - T_W
            else:
                wait_time = records[i][1] - records[i - 1][2]
                W_sum += wait_time
                if wait_time - T_W > 0:
                    w_thanT_sum += wait_time - T_W
    fitn = max_F + W_sum + w_thanT_sum * gamma
    return fitn, max_F, W_sum, w_thanT_sum


def get_max_makespan(service_status):
    F = []
    for i in range(1, len(service_status)):
        recs = service_status[i]
        recs.sort(key=lambda a: a[2])
        F.append(recs[-1][2])
    return max(F)


def get_project_bound(p, m):
    """
    :param p 客户编号 [1, n]
    :param m 检查项目的数目
    """
    start = (p - 1) * m + 1
    end = p * m
    return [i for i in range(start, end + 1)]


def gen_patient_proj_idx_seq(chromosome, bound):
    """
    :param chromosome: 染色体
    :param bound: 顾客项目编码范围
    :return: 两个list，一个是顾客项目编码在染色体中的index，一个是顾客项目编码值
    """
    idx = [i for i in range(len(chromosome.genes)) if chromosome.genes[i] in bound]
    seq = [val for val in chromosome.genes if val in bound]
    return idx, seq


def change_patient(chromosome, own_change_index, other_proj):
    """
    换某一顾客的全部项目
    :param chromosome: 染色体备份-浅拷贝就够用
    :param own_change_index: 自己要更换顾客项目的各个位置index
    :param other_proj: deque结构，要更换的项目，要求保持之前的顺序
    :return: 新的染色体
    """
    for i in range(len(chromosome.genes)):
        if i in own_change_index:
            chromosome.genes[i] = other_proj.pop(0)  # 从左边出队
    return chromosome


class SGA:
    def __init__(self, people_num, service_times, cross_rate, mutation_rate, pop_size):
        self.people_num = people_num
        self.service_times = service_times
        self.project_num = len(service_times) - 1
        self.pc = cross_rate
        self.pm = mutation_rate
        self.pop_size = pop_size
        self.pop = self.__init_pop()

    def __init_pop(self):
        pop = []
        for i in range(self.pop_size):
            pop.append(DNA(self.people_num, self.service_times))
        return pop

    def translate(self):
        """解码,之后可以获得种群的fitness之和"""
        for dna in self.pop:
            self.translate_dna(dna)

    def translate_dna(self, dna):
        """单条dna解码，保存状态并计算适应度"""
        if dna.patient_records.not_full() is False:
            raise ValueError("存在重复记录!")
        # 创建记录
        for gene in dna:
            self.create_records(gene, dna)
        # 计算指标
        dna.fitness = compute_fitness(dna)

    def check_record(self, records, record):
        """检查相同记录"""
        for rec in records:
            if rec[1] == record[1]:
                return True

    def filter_checked(self, records, times):
        """过滤可能造成相同记录的时间段"""
        if len(times) == 0:
            return []
        filtered_times = [ti for ti in times if self.check_record(records, ti) is False]
        return filtered_times

    def create_records(self, gene, dna):
        """解码后为客户和服务台分别创建一条记录"""
        pidx, sidx = self.get2index(gene)
        consume_time = self.service_times[sidx]

        # 插入记录
        p_status = dna.patient_records  # 当前dna的客户表格
        ser_status = dna.service_records  # 当前dna的服务台表格

        last_end_time = p_status.get_last_project_end(pidx)  # 上一项目结束时间
        idles = ser_status.get_idle_times(sidx)  # 当前服务台sidx的空闲区间
        if len(idles) == 1 and idles[0][0] == 0 and idles[0][1] == TIME_MAX:  # 科室没有人被分配
            cur_start_time = max(last_end_time, 0)
            cur_end_time = cur_start_time + consume_time  # 得到结束时间
            p_status.add_record(pidx, [sidx, cur_start_time, cur_end_time])
            try:
                ser_status.add_record(sidx, [pidx, cur_start_time, cur_end_time])
            except ValueError:
                print("第一条数据有错!")

        else:
            # 空闲时间符合插入记录的

            enough = [time for time in idles if time[1] - time[0] > consume_time]
            enough = self.filter_checked(ser_status.store[sidx], enough)
            # 根据起始时间排序
            enough.sort(key=lambda e: e[0])
            if len(enough) > 0:  # 有足够时间
                idle = enough[0]
                cur_start_time = max(last_end_time, idle[0])
                cur_end_time = cur_start_time + consume_time
                p_status.add_record(pidx, [sidx, cur_start_time, cur_end_time])
                ser_status.add_record(sidx, [pidx, cur_start_time, cur_end_time])

            else:  # 没有足够时间, 设置有超时的
                records = ser_status.store[sidx]
                records.sort(key=lambda a: a[2])
                cur_start_time = max(last_end_time, records[-1][2])
                cur_end_time = cur_start_time + consume_time
                p_status.add_record(pidx, [sidx, cur_start_time, cur_end_time])
                ser_status.add_record(sidx, [pidx, cur_start_time, cur_end_time])

    def get_fitness(self):
        """获取整个种群的fitness"""
        all_fit = 0
        fits = []
        for dna in self.pop:
            fits.append(dna.fitness)
            all_fit += dna.fitness
        return all_fit, np.array(fits)

    def get2index(self, item):
        """得到编码解码后客户下标和服务台下标"""
        div = item // self.project_num
        mod = item % self.project_num
        if mod == 0:
            pidx = div
            sidx = self.project_num
        else:
            pidx = div + 1
            sidx = mod
        return pidx, sidx

    def select(self):
        """随机选择交叉Dna"""
        fit_sum, fits = self.get_fitness()
        pop_copy = deepcopy(self.pop)
        # pop_copy.sort(key=lambda a:a.fitness)
        # pop_copy = pop_copy[:16]
        # select_pop = []
        # i = 0
        # while i < self.pop_size:
        #     idx = np.random.randint(len(pop_copy))
        #     select_pop.append(pop_copy[idx])
        #     i += 1
        # assert len(select_pop) == self.pop_size
        # best_dna = deepcopy(self.pop[np.argmin(fits)])
        # return select_pop, best_dna
        select_pop = []
        for i in range(self.pop_size):
            idx = np.random.randint(len(pop_copy))
            select_pop.append(pop_copy[idx])
        assert len(select_pop) == self.pop_size
        best_dna = deepcopy(self.pop[np.argmin(fits)])
        return select_pop, best_dna

    def crossover(self, P1, P2):
        """交叉"""
        if np.random.rand() < self.pc:
            cross_idx = np.random.randint(1, self.people_num + 1)  # 随机获取一个顾客的编号
            project_num = self.project_num
            change_bound = get_project_bound(cross_idx, project_num)  # 解码顾客项目编号范围
            # 获取两条染色体上的顾客项目idx和seq
            idx1, seq1 = gen_patient_proj_idx_seq(P1, change_bound)
            idx2, seq2 = gen_patient_proj_idx_seq(P2, change_bound)
            # 根据一定概率pc交叉产生子代C1， C2
            C1 = change_patient(deepcopy(P1), idx1, seq2)
            self.mutation(C1)
            C2 = change_patient(deepcopy(P2), idx2, seq1)
            self.mutation(C2)

            # self.translate_dna(P1)
            C1.clear_records()
            C2.clear_records()
            self.translate_dna(C1)
            self.translate_dna(C2)
            fits = [P1.fitness, C1.fitness, C2.fitness]
            idx = np.argmin(fits)
            if idx == 0:
                return P1
            elif idx == 1:
                return C1
            else:
                return C2
        else:
            return P1

    @staticmethod
    def mutation(chromosome):
        """
        基因突变，任意两个基因交换顺序
        :param chromosome: 染色体
        """
        bound = len(chromosome.genes)
        idx1, idx2 = np.random.choice([i for i in range(0, bound)], size=2, replace=False)
        chromosome.genes[idx1], chromosome.genes[idx2] = chromosome[idx2], chromosome[idx1]

    def evolve(self):
        """将交叉变异后的优良染色体引入种群"""
        select_pop, best_dna = self.select()
        new_pop = []
        for parent in select_pop:  # for every parent
            child = self.crossover(parent, best_dna)
            new_pop.append(child)
        self.pop = new_pop


def save_dna(gen_str, dna, metric_str, path):
    dna_seq = dna.__repr__()
    lines = [gen_str, metric_str, dna_seq]
    with open(path, "a+", encoding='utf-8') as f:
        f.write('\n'.join(lines))


def run(algo, iter_num):
    best_fits = []
    temp_dna = DNA(TOTAL_PEOPLE_NUM, FIXED_SERVICE_TIMES)
    temp_dna.fitness = 1000000
    temp_fits = 1000000
    for generation in range(iter_num):
        # 解码种群
        algo.translate()
        # 计算适应度值
        _, fits = algo.get_fitness()
        # 拿到适应度值最小的dna
        best_idx = np.argmin(fits)
        # 输入或者记录到list——fitness
        gen_str = '\nGen: {} | best fit: {}'.format(generation, fits[best_idx])
        print(gen_str)
        best_fits.append(fits[best_idx])
        d = algo.pop[best_idx]
        if d.fitness < temp_fits:
            temp_fits = d.fitness
            temp_dna = d
            temp_gen = gen_str
            metric = compute_metric(temp_dna)
            metric_str = 'fitness: {}, maxF: {}, W_sum: {}, w_thanT_sum: {}'.format(metric[0], metric[1], metric[2],
                                                                                    metric[3])
            save_dna(temp_gen, temp_dna, metric_str, './result2.txt')

        # 优良染色体加入种群
        algo.evolve()
        for dna in algo.pop:
            dna.clear_records()
    return best_fits



if __name__ == '__main__':
    POPULATION_SIZE = 150  # 种群大小, 即解的个数
    N_GENERATIONS = 550  # 迭代次数
    CROSS_RATE = 0.8  # 交叉概率
    MUTATE_RATE = 1.0  # 变异概率
    FIXED_SERVICE_TIMES = [
        -1,  # 0无
        3,  # 1体质测试
        3,  # 2内科
        4,  # 3外科
        2,  # 4眼耳口鼻科
        3,  # 5验血
        2,  # 6心电图
        5,  # 7X光
        6,  # 8B超
    ]  # 共28分钟
    TOTAL_PEOPLE_NUM = 10

    TEST_PEOPLE_NUM = 10
    TEST_POP_SIZE = 1
    TEST_SERVER_TIMES = [
        -1,
        2,
        2,
        6
    ]
    params = {
        'people_num': TOTAL_PEOPLE_NUM,
        'service_times': FIXED_SERVICE_TIMES,
        'cross_rate': CROSS_RATE,
        'mutation_rate': MUTATE_RATE,
        'pop_size': POPULATION_SIZE,
    }

    test_params = {
        'people_num': TEST_PEOPLE_NUM,
        'service_times': TEST_SERVER_TIMES,
        'cross_rate': CROSS_RATE,
        'mutation_rate': MUTATE_RATE,
        'pop_size': TEST_POP_SIZE,
    }

    sga = SGA(**params)
    start = time.time()
    data = run(sga, N_GENERATIONS)
    end = time.time()
    print("运行时间: {}s".format(end - start))
    plt.plot(data)
    plt.show()
