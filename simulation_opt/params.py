cost_time_lookup = [
    3,  # 0体质测试
    3,  # 1内科
    4,  # 2外科
    2,  # 3眼耳口鼻科
    3,  # 4验血
    2,  # 5心电图
    5,  # 6X光
    6,  # 7B超
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
total_people = 30
project_num = len(cost_time_lookup)
T_W = 15  # 等待阈值
GAMMA1 = 1  # 惩罚系数1，用于实验1
GAMMA2 = 10  # 惩罚系数2，用于超阈值等待
POP_SIZE = 10  # 种群大小
GROUP = 30  # 选择权重，百分之百
CROSS_RATE = 0.8
MUTATE_RATE = 1.0
N_GENERATIONS = 500