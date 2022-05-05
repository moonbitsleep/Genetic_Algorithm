import numpy as np
"""
======================================
    自由到达的顾客在系统开放前30min遵循参数2lam的泊松分布，
    随后的150分钟遵循参数为lam的泊松分布
    启发式1: 选择人数最少的队列
    启发式2: 贪心最短
=======================================    
"""


if __name__ == '__main__':
    print(np.random.poisson(lam=114, size=60))