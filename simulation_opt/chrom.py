from params import *
import random
import copy

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

    def scan_service_table(self):
        """解码后局部搜索使用"""
        self.compute_fitness()
        for i in range(len(self.service_table)):
            queues = self.service_table[i]
            for j in range(len(queues)):
                queue_records = queues[j]
                for k in range(len(queue_records) - 2):
                    record1 = queue_records[k]
                    record2 = queue_records[k + 1]
                    if record2[1] - record1[2] >= 15:
                        return i, record2[0], queue_records[k+2][0]
        proj_idx = random.randint(0, len(self.service_table) - 1)
        que_idx = random.randint(0, len(self.service_table[proj_idx]) - 1)
        rec_idx1 = random.randint(0, len(self.service_table[proj_idx][que_idx]) - 2)
        rec_idx2 = rec_idx1 + 1
        return proj_idx, self.service_table[proj_idx][que_idx][rec_idx1][0], self.service_table[proj_idx][que_idx][rec_idx2][0]

    def explore_search1(self):
        proj, p1, p2 = self.scan_service_table()
        gene1 = self.deserialization(proj, p1)
        gene2 = self.deserialization(proj, p2)
        seq_copy = copy.deepcopy(self.sequence)
        idx1 = seq_copy.index(gene1)
        idx2 = seq_copy.index(gene2)
        seq_copy[idx1], seq_copy[idx2] = seq_copy[idx2], seq_copy[idx1]
        return Chromosome(seq_copy)

    @staticmethod
    def deserialization(proj_idx, pat_idx):
        """根据项目index和顾客index得到gene"""
        gene = proj_idx * project_num + pat_idx + 1
        return gene

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
        找到某科室最早空闲的服务台
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
