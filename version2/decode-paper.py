import copy

cost_time_lookup = [
    3,  # 0体质测试   4
    3,  # 1内科      4
    4,  # 2外科      3
    # 2,  # 3眼耳口鼻科 5
    # 3,  # 4验血      4
    # 2,  # 5心电图    5
    # 5,  # 6X光      3
    6,  # 7B超      2
]  # 共28分钟


total_people = 12
project_num = len(cost_time_lookup)
T_W = 15
GAMMA = 5

def translate_operation(opt):
    """解码操作 （patient_index, project_index）"""
    opt = opt - 1
    return opt // project_num, opt % project_num

def translate_seq(seq):
    """解码整个序列"""
    res = []
    for o in seq:
        res.append(translate_operation(o))
    print(res)

def get_people_last_end(people_table):
    if len(people_table) == 0:
        return 0
    else:
        people_table.sort(key=lambda e: e[2])
        return people_table[-1][2]


def get_project_last_end(project_table):
    if len(project_table) == 0:
        return 0
    else:
        project_table.sort(key=lambda e: e[2])
        return [record[-1] for record in project_table]

def get_middle_start_time(end_time_list, cost_time, people_start):
    for i in range(1, len(end_time_list)):
        if end_time_list[i] - end_time_list[i - 1] >= 2 * cost_time:
            if end_time_list[i] - people_start >= 2 * cost_time:
                return max(people_start, end_time_list[i - 1])
    return -1


def build_table(seq):
    global cost_time_lookup
    global total_people
    global project_num
    people_records = [[] for _ in range(total_people)]
    project_records = [[] for _ in range(project_num)]
    for operation in seq:
        people_index, project_index = translate_operation(operation)
        cost_time = cost_time_lookup[project_index]

        people = people_records[people_index]  # 客户的所有记录
        project = project_records[project_index]  # 项目的所有记录
        # 创建记录
        people_last_end_time = get_people_last_end(people)  # 人有空的最早时间
        project_last_ends = get_project_last_end(project)

        # 科室第一条记录
        if type(project_last_ends) == int:
            start_time = people_last_end_time
            end_time = start_time + cost_time
            people.append([project_index, start_time, end_time])
            project.append([people_index, start_time, end_time])

        if type(project_last_ends) == list:

            tmp_start = get_middle_start_time(project_last_ends, cost_time, people_last_end_time)
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


def compute_metrics(seq):
    """根据dna计算指标"""
    global GAMMA
    global T_W
    people_records, project_records = build_table(seq)
    W_sum = 0
    w_thanT_sum = 0
    tmp_lates = []
    for records in people_records:
        records.sort(key=lambda e: e[2])
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
    fitness = maxF + W_sum + GAMMA * w_thanT_sum
    print("****************")
    print("dna seq:", seq)
    print("people's records: ")
    for p in people_records:
        print(p)
    print("service records: ")
    for s in project_records:
        s.sort(key=lambda e: e[2])
        print(s)
    print("fitness:", fitness)
    print("makespan:", maxF)
    print("total_wait:", W_sum)
    print("greater_than_threshold:", w_thanT_sum)
    print("****************")
    return fitness, maxF, W_sum, w_thanT_sum


if __name__ == '__main__':
    dna = [46, 16, 34, 7, 3, 47, 32, 9, 22, 33, 23, 8, 1, 45, 18, 27, 21, 36, 4, 35, 40, 14, 25, 26, 15, 11, 19, 24, 48, 42, 31, 38, 44, 43, 37, 39, 28, 30, 12, 5, 29, 13, 41, 10, 17, 20, 2, 6]

    compute_metrics(dna)

    better_dna = [43, 5, 9, 13, 40, 33, 12, 7, 41, 8, 30, 44, 29, 2, 17, 46, 26, 31, 25, 39, 19, 22, 6, 47, 10, 42, 1, 15, 32, 38, 37, 20, 16, 14, 3, 21, 11, 4, 27, 45, 48, 34, 28, 35, 23, 36, 18, 24]
    compute_metrics(better_dna)

