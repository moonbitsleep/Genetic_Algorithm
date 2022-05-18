import numpy as np


def OCBA(design_mean, design_var, design_used, addiction_budget):
    """
    The first round is not considered in this function, you need to set initional budget
    for each design and simulate in the calling funcion. When you get the static information
    after simulation, you can all this function to calculate budget needed for each design in next round.

    Points:
    1. Must ensure that the addiction_budget is all allocated, no more than or no remaining
    2. Dont consider participation for the design with var 0
    3. If there are multiple optimal designs (which have maximum mean),
        select the first appeared optimal design to be
        the best design and select the first appeared design which dont has maximum mean
        to be base design for calculating radio among all designs
    :param design_mean: Mean  of each design's performance,
        must be numpy.ndarray for data type and (N,) for data shape
    :param design_var: Variance of each design's performance,
        must be numpy.ndarray for data type and (N,) for data shape
    :param design_used: Computing budget already used by each design,
        must be numpy.ndarray for data type and (N,) for data shape
    :param addiction_budget: Newly added computing budget,
        must be numpy.ndarray or list for data type and (N,) for data shape
    :return: Computing budget allocated to each design with addiction_budget given in next round
    """

    assert addiction_budget > 0, "wrong addiction_budget"
    assert len(design_mean) >= 2, "amount of designs must greater than 1"
    assert isinstance(design_mean, np.ndarray) and isinstance(design_var, np.ndarray) and isinstance(design_used, np.ndarray), \
        "type of design_mean,design_var,design_used must be numpy.ndarray"

    def Normalization(x):
        sum_ = sum(x)
        x /= sum_
        return x


    design_num = len(design_mean)
    no_zero_var_ids = [i for i in range(design_num) if design_var[i] != 0.]
    best_ids = list(np.where(design_mean[no_zero_var_ids] == max(design_mean))[0])
    other_ids = [i for i in no_zero_var_ids if i not in best_ids]
    var_highest_id = np.argmax(design_var)

    # If the var of the best design is 0, then OCBA ends and returns False
    # and you can randomly selects an optimal action in the calling function
    if not best_ids or len(no_zero_var_ids) == 1:
        return False, 0

    radio = np.zeros(len(design_mean))
    best_mean = design_mean[best_ids[0]]

    # compute radio of every design
    # Discuss in different situations:
    # whether all designs are best design, if all, calculate the radio according to the variance of each design
    # (someone may randomly choose one design to have radio 1 in order to allocate all addiction budget to it)
    # if not all, normal calculate the radio according to the mean and variance information
    if len(best_ids) == len(no_zero_var_ids):  # all designs are best
        radio[best_ids[0]] = 1.
        radio[best_ids] = design_var[best_ids] / design_var[best_ids[0]]
    else:
        base_id = other_ids[0]
        base_mean = design_mean[base_id]
        base_var = design_var[base_id]
        radio[base_id] = 1.
        temp = 1. / design_var[base_id]
        for id in other_ids[1:]:
            radio[id] = (design_var[id] * np.square(best_mean - base_mean)) \
                / (base_var * np.square(best_mean - design_mean[id]))
            temp += np.square(radio[id]) / design_var[id]
        radio[best_ids] = np.sqrt(design_var[best_ids] * temp)
    radio = Normalization(radio)

    # Cycically calculate the computing budget that each design should have used
    active = np.array([False] * len(design_mean))
    active[no_zero_var_ids] = True
    should_used = design_used.copy()
    any_exceed = True
    while(any_exceed):
        temp_budget = sum(design_used[active]) + addiction_budget
        should_used[active] = np.floor(temp_budget * radio)[active]
        exceed = should_used < design_used
        any_exceed = any(exceed)
        if any_exceed:
            should_used[exceed] = design_used[exceed]
            active[exceed] = False
            radio[~active] = 0
            radio = Normalization(radio)

    # Allocating unused computing budget to the design with the highest variance
    should_used[var_highest_id] += temp_budget - sum(should_used[active])
    assert sum(should_used - design_used) == addiction_budget, "wrong"
    
    return True, should_used - design_used


if __name__ == '__main__':
    mean = np.arange(1, 501)
    var = np.arange(1, 501)
    used = np.ones(500) * 30
    budget = 66000
    flag, ocba = OCBA(mean, var, used, budget)
    print(int(ocba[-1]))
