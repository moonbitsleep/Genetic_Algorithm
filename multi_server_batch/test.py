from collections import deque
from operator import attrgetter

import numpy as np


class A:
    def __init__(self, idx, age):
        self.id = idx
        self.age = age

    def __repr__(self):
        return str(self.age)

if __name__ == '__main__':
    # a = deque([A(0, 6), A(1, 1), A(2, 2), A(3, 19), A(4, 9), A(5, 19)])
    # b = sorted(a, key=attrgetter("age"))
    # print(b)
    # b.insert(2, A(6, 5))
    # print(b)
    a = np.random.choice([i for i in range(64)], size=64, replace=False)
    b = deque(a)
    print(b)
    if 9 in b:
        print("true")
    else:
        print("None")
