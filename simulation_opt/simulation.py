import math
import random
import time
from params import *
import numpy as np
import matplotlib.pyplot as plt
from chrom import Chromosome
from population import Population
"""
======================================
    根据期望时间生成一个随机服务时间的矩阵
    行数 = 人数
    列数 = 项目数
=======================================    
"""

def generate_time_table(people_total, project_total):
    res = np.zeros(shape=(project_total, people_total))
    for i in range(len(cost_time_lookup)):
        res[i, :] = np.random.exponential(cost_time_lookup[i], [1, people_total])

    res = res.flatten().astype(int)
    for i in range(len(res)):
        if res[i] == 0:
            res[i] = 1
    res = res.reshape(project_total, people_total)
    return res.T

def SGA():
    population = Population(POP_SIZE)


if __name__ == '__main__':
    pass

