import numpy as np

cost_time_lookup = [
    3,
    3,
    4,
    2,
    3,
    2,
    5,
    6,
]

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



if __name__ == '__main__':
    pass