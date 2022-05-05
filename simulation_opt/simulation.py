import math
import random
import numpy as np
import matplotlib.pyplot as plt

"""
======================================
    根据期望时间生成一个随机服务时间的矩阵
    行数 = 人数
    列数 = 项目数
=======================================    
"""
expect_service_times = [
    3,  # 0体质测试
    3,  # 1内科
    4,  # 2外科
    2,  # 3眼耳口鼻科
    3,  # 4验血
    2,  # 5心电图
    5,  # 6X光
    6,  # 7B超
]


def expntl(L):
    """
    negative exponential distribution
    return a double random number, L is the mean value
    """
    u = random.random()
    return -L * math.log(u)


def gen_service_time_matrix():
    pass


if __name__ == '__main__':
   total = 60
   serve_times = [math.ceil(expntl(6)) for _ in range(total)]
   print(serve_times)
   print(np.mean(serve_times))
   plt.hist(serve_times)
   plt.show()