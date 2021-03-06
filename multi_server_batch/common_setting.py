"""
====================================
            exp params
====================================
"""
import math

cost_time_lookup = [
    3,
    3,
    4,
    2,
    3,
    2,
    5,
    6,
]  # 28 mins
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

TOTAL_PEOPLE = 84
PROJECT_NUM = len(cost_time_lookup)
END_TIME = 240
BATCH_INTERVAL_TIME = 30
BATCH_COUNT = (END_TIME - 60) // BATCH_INTERVAL_TIME  + 1 # 7
BATCH_PEOPLE = math.ceil(TOTAL_PEOPLE / BATCH_COUNT) # n
T_W = 15  # wait threshold
GAMMA1 = 1  # penalty coefficient for normal
GAMMA2 = 10  # penalty coefficient for wait threshold

"""
============================================
                genetic algorithm
============================================
"""
POP_SIZE = 100
TOTAL_POP = 10  # one for save others for evolve
CROSS_RATE = 0.8
MUTATE_RATE = 0.4
GEN_MAX = 30  # The continuous GEN_MAX generation retains the previous optimal value

if __name__ == '__main__':
    print(TOTAL_RESOURCE)