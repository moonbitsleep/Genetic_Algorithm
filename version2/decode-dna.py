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
    """得到项目表中每个人的结束时间"""
    if len(project_table) == 0:
        return 0
    else:
        project_table.sort(key=lambda e: e[2])
        return project_table[-1][2]

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

        start_time = max(people_last_end_time, project_last_ends)
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
    return fitness, maxF, W_sum, w_thanT_sum


if __name__ == '__main__':
    d1 = [29, 43, 13, 23, 5, 28, 46, 35, 11, 14, 39, 21, 36, 24, 27, 34, 47, 38, 33, 41, 8, 7, 37, 12, 44, 9, 30, 19, 42, 18, 45, 40, 25, 1, 32, 26, 31, 16, 17, 2, 22, 48, 6, 15, 20, 4, 3, 10]
    d2 = [17, 36, 31, 24, 1, 14, 3, 16, 6, 30, 38, 43, 29, 45, 26, 7, 10, 35, 20, 34, 42, 39, 5, 8, 32, 2, 25, 27, 21, 9, 15, 19, 23, 28, 37, 11, 13, 4, 46, 18, 33, 40, 41, 47, 44, 12, 48, 22]

    print(compute_metrics(d1))
    d1_ = copy.deepcopy(d1)

    print(compute_metrics(d1_))