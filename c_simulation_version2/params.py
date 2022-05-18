import math

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
TOTAL_RESOURCE = sum(resource_look_up)

TOTAL_PEOPLE = 60
PROJECT_NUM = len(resource_look_up)
END_TIME = 240
BATCH_INTERVAL_TIME = 30
BATCH_COUNT = (END_TIME - 60) // BATCH_INTERVAL_TIME  + 1 # 7
BATCH_PEOPLE = math.ceil(TOTAL_PEOPLE / BATCH_COUNT) # n
T_W = 15  # wait threshold
GAMMA1 = 1  # penalty coefficient for normal
GAMMA2 = 10  # penalty coefficient for wait threshold


POP_SIZE = 100
TOTAL_POP = 10  # one for save others for evolve
CROSS_RATE = 0.8
MUTATE_RATE = 0.4
GEN_MAX = 30  # The continuous GEN_MAX generation retains the previous optimal value