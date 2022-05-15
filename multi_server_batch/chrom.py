from common_setting import *
from collections import deque

class Chromosome:
    def __init__(self, seq):
        self.seq = seq
        self.fitness = None
        self.makespan = None
        self.total_wait = None
        self.greater_than_threshold = None
        self.service_idle = None
        self.patients_table = [
            [deque() for _ in range(BATCH_PEOPLE)] for _ in range(BATCH_COUNT)
        ]
        self.services_table = [
            [deque() for _ in range(resource_look_up[sidx])] for sidx in range(PROJECT_NUM)
        ]

    @staticmethod
    def translate_gene(gene):
        gene = gene - 1
        return gene // PROJECT_NUM, gene % PROJECT_NUM

    @staticmethod
    def get_people_last_end(people_records, batch_idx):
        if len(people_records) == 0:
            return batch_idx * BATCH_INTERVAL_TIME
        else:
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
                if queue[-1][2] < targets[1]:
                    targets[0] = i
                    targets[1] = queue[-1][2]
        return targets[0], targets[1]

    @staticmethod
    def add_queue_records(queue_records, record):
        flag = len(queue_records)
        for i in range(len(queue_records) - 1):
            if queue_records[i][2] <= record[1] and queue_records[i+1][1] >= record[2]:
                flag = i
        if flag == len(queue_records):
            queue_records.append(record)
        else:
            queue_records.insert(flag+1, record)

    def translate(self):
        for batch in range(len(self.seq)):
            sequence = self.seq[batch]
            for gene in sequence:
                peo_idx, proj_idx = self.translate_gene(gene)
                cost_time = cost_time_lookup[proj_idx]
                peo_records = self.patients_table[batch][peo_idx]
                queues_records = self.services_table[proj_idx]

                peo_early_time = self.get_people_last_end(peo_records, batch)
                que_idx, proj_early_time = self.get_early_queue(queues_records)
                start_time = max(peo_early_time, proj_early_time)
                end_time = start_time + cost_time

                #add
                peo_records.append([proj_idx, que_idx, start_time, end_time])
                self.add_queue_records(queues_records[que_idx], [peo_idx, start_time, end_time])

    def clean_table_cache(self):
        self.patients_table = [
            [deque() for _ in range(BATCH_PEOPLE)] for _ in range(BATCH_COUNT)
        ]
        self.services_table = [
            [deque() for _ in range(resource_look_up[sidx])] for sidx in range(PROJECT_NUM)
        ]

    def compute_metric(self):
        if self.fitness is not None:
            self.clean_table_cache()
        self.translate()
        W_sum = 0
        w_thanT_sum = 0
        tmp_end = deque()
        for batch in range(len(self.patients_table)):
            batch_records = self.patients_table[batch]
            for records in batch_records:
                assert len(records) == PROJECT_NUM
                tmp_end.append(records[-1][3])
                for i in range(len(records)):
                    if i == 0:
                        wait_time = records[i][2] - batch * BATCH_INTERVAL_TIME
                        assert wait_time >= 0
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

        idle_sum = 0
        for proj_idx in range(len(self.services_table)):
            queues_records = self.services_table[proj_idx]
            for records in queues_records:
                for i in range(len(records)):
                    if i == 0:
                        idle_time = records[i][1] - 0
                        if idle_time > 0:
                            idle_sum += idle_time
                    else:
                        idle_time = records[i][1] - records[i-1][2]
                        if idle_time > 0:
                            idle_sum += idle_time

        self.makespan = max(tmp_end)
        self.total_wait = W_sum
        self.greater_than_threshold = w_thanT_sum
        self.service_idle = idle_sum
        self.fitness = GAMMA1 * (self.total_wait + self.service_idle) + GAMMA2 * (self.greater_than_threshold + self.makespan - END_TIME)


    def print_metric(self):
        print("fitness: ", self.fitness)
        print("over work time: ", self.makespan - END_TIME)
        print("avg_wait: ", self.total_wait / TOTAL_PEOPLE)
        print("greater_than_threshold: ", self.greater_than_threshold)
        print("avg service_idle: ", self.service_idle / TOTAL_RESOURCE)

    def print_service_records(self):
        for service in self.services_table:
            for records in service:
                print(records)

    def print_people_records(self):
        for batch in self.patients_table:
            for records in batch:
                print(records)


if __name__ == '__main__':
    c = Chromosome([1, 2, 3])


