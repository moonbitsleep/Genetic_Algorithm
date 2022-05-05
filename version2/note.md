# 前置

## 概述

遗传算法在全局搜索中具有较强的性能，但在求解约束条件复杂且解空间较大的问题时，存在
局部搜索能力()较差、易早熟收敛且消耗时间长等不足。   
改进：62, 63
- 算法编码解码方法
- 种群初始化
- 遗传算子的设计和算法逻辑设计
  多种群遗传算法: 64, 65, 66  
  

离散随机优化问题：(1)解空间抽样和搜索;(2)评估样本解的性能值.
仿真优化方法67
序优化和最佳计算预算分配方法68
69-72

## 论文参考

62 王雷 基于改进遗传算法的柔性作业车间调度
63 Chien C F, An evolutionary approach to rehabilitation patient scheduling: A case Study.
64 蔡良伟 作业车间调度问题的多种群遗传算法
65 用带蚁群搜索的多种群遗传算法求解作业车间调度问题
66 多目标模糊柔性车间调度中的多种群遗传算法
67 仿真优化理论和方法综述
69 An ordinal optimization-based evolution strategy to schedule complex make-to-order products
70 Evolutionary algorithm for stochastic job shop scheduling with random processing time.
71 Optimal Computing Budget Allocation for Ordinal Optimization in Solving Stochastic Job Shop Scheduling Problems
72 A Novel Stochastic Simulation Optimization Method in Solving Job Shop Scheduling Problem Under Processing Time Variability
# 确定服务时长+单种群+局部搜索
## 问题描述
确定每个项目上的顾客顺序，确定每个顾客的项目顺序
调度轨迹[(pid, sid), (pid, sid)]
实验设置：
1. 顾客在体检机构开放的0时刻全部到达
1. 顾客全部参加重点科室的检查，每种科室对每位顾客的服务时间相同。
1. 每种科室资源仅有一个服务台。
1. 各科室检查顺序无约束。

实验指标metric
$$
\min \left(\sum_{i=1}^{n} W_{i}+\max _{i=1, \ldots, n}\left(F_{i}\right)+\gamma \times \sum_{i=1}^{n} \sum_{j=1}^{m}\left(w_{i j}-T^{W}\right)^{+}\right)
$$

## 混合遗传算法

步骤1):随机产生初始解，建立初始种群，种群规模为Z，进行**<span style="color:red;">解码和适应度计算</span>**。

步骤2):对初始种群中的每个解应用局部搜索程序，使其改进为局部最优(涉及**<span style="color:red;">解码和适应度计算</span>**)。

步骤3):使用交叉算子和突变算子对当前种群中的解进行重组，产生后代，形成新的种群。进行**<span style="color:red;">解码和适应度计算</span>**。

步骤4):对重组后的种群中的每个解应用局部搜索程序，使其改进为局部最优(涉及**<span style="color:red;">解码和适应度计算</span>**)。

步骤5):重复3-4，直到解一定代数内最优值没有发生明显变化改进或到达最大迭代次数，输出当前最优解。

### 染色体编码和解码

```python
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

total_people = 20
project_num = len(cost_time_lookup)

def translate_operation(opt):
    """解码操作"""
    opt = opt - 1
    return opt // project_num, opt % project_num


class Chromosome:
    def __init__(self, sequence):
        self.sequence = sequence  # list DNA序列
        self.fitness = None
        self.makespan = None
        self.maxF = None
        self.total_wait = None
        self.greater_than_threshold = None

    def translate(self):
        global cost_time_lookup
        global total_people
        global project_num
        people_records = [[] for _ in range(total_people)]
        project_records = [[] for _ in range(project_num)]
        for operation in self.sequence:
            people_index, project_index = translate_operation(operation)
            cost_time = cost_time_lookup[project_index]

            people = people_records[people_index] # 客户的所有记录
            project = project_records[project_index] # 项目的所有记录
            # 创建记录
            people_last_end_time = self.get_people_last_end(people) # 人有空的最早时间
            project_last_ends = self.get_project_last_end(project)

            # 科室第一条记录
            if type(project_last_ends) == int:
                start_time = people_last_end_time
                end_time = start_time + cost_time
                people.append([project_index, start_time, end_time])
                project.append([people_index, start_time, end_time])

            if type(project_last_ends) == list:
            # 科室中间插入
                pass
            # 科室末尾插入



    @staticmethod
    def get_people_last_end(people_table):
        if len(people_table) == 0:
            return 0
        else:
            people_table.sort(key=lambda e: e[2])
            return people_table[-1][2]

    @staticmethod
    def get_project_last_end(project_table):
        """得到项目表中每个人的结束时间"""
        if len(project_table) == 0:
            return 0
        else:
            project_table.sort(key=lambda e: e[2])
            return [record[-1] for record in project_table]

    @staticmethod
    def get_middle_start_time(end_time_list, cost_time, people_start):
        """
        如果二者间时间间隔大于2倍的检查时间，表明有足够的时间
        如果后者 - 客户开始时间仍大于2倍的检查时间(要减去之前分配的一次检查)，则表明可以插入
        否则，只能在末尾插入
        """
        for i in range(1, len(end_time_list)):
            if end_time_list[i] - end_time_list[i - 1] >= 2 * cost_time:
                if end_time_list[i] - people_start >= 2 * cost_time:
                    return max(people_start, end_time_list[i - 1])
        return -1

if __name__ == '__main__':
    seq = [8, 9, 24, 26, 36]
    c = Chromosome(seq)
    c.translate()

    pro_table = [[2, 0, 2], [10, 2, 4], [6, 12, 14], [1, 14, 16], [3, 24, 26], [5, 42, 44], [4, 44, 46], [9, 54, 56], [8, 56, 58], [7, 60, 62]]
    ta = []
    a = c.get_project_last_end(pro_table)
    start = c.get_middle_start_time(a, 2, 50)
    print(a)
    print(start)
    [2, 4, 14, 16, 26, 44, 46, 56, 58, 62]
	50
    
```





### 局部搜索

```markdown
Gen: 322, best dna seq: [169, 21, 127, 212, 201, 43, 29, 210, 224, 116, 56, 9, 36, 96, 120, 209, 151, 223, 82, 94, 66, 154, 236, 225, 125, 23, 86, 129, 63, 132, 220, 108, 34, 59, 124, 174, 85, 55, 211, 12, 218, 222, 165, 235, 5, 149, 49, 39, 117, 73, 6, 30, 20, 163, 214, 44, 202, 193, 100, 78, 208, 105, 87, 11, 53, 122, 16, 138, 180, 166, 194, 187, 4, 48, 134, 22, 190, 98, 95, 28, 186, 145, 109, 69, 18, 92, 130, 133, 52, 157, 81, 228, 200, 62, 238, 24, 68, 51, 102, 195, 61, 70, 50, 114, 213, 121, 26, 60, 234, 10, 65, 35, 42, 215, 128, 140, 58, 171, 217, 110, 126, 197, 77, 19, 156, 182, 175, 27, 17, 204, 173, 38, 147, 113, 1, 91, 74, 219, 101, 90, 150, 207, 123, 15, 57, 137, 106, 206, 88, 170, 45, 40, 93, 198, 47, 139, 199, 146, 118, 185, 41, 37, 46, 161, 76, 13, 216, 131, 83, 176, 230, 64, 7, 103, 205, 2, 152, 221, 25, 99, 136, 14, 233, 158, 141, 189, 119, 203, 3, 177, 89, 142, 148, 111, 226, 167, 67, 153, 31, 188, 143, 135, 104, 75, 192, 196, 168, 191, 97, 33, 54, 8, 229, 71, 144, 32, 79, 239, 112, 181, 227, 107, 240, 178, 115, 231, 172, 162, 237, 84, 155, 179, 72, 232, 164, 184, 160, 80, 159, 183]
best fitness: 5658
makespan: 180, total_wait: 2533, T_W_wait: 589
Gen: 323, best dna seq: [169, 21, 127, 212, 201, 43, 29, 210, 224, 116, 56, 9, 36, 96, 120, 209, 151, 223, 82, 94, 66, 154, 236, 225, 125, 23, 86, 129, 63, 132, 220, 108, 34, 59, 124, 174, 85, 55, 211, 12, 218, 222, 165, 235, 5, 149, 49, 39, 117, 73, 6, 30, 20, 163, 214, 44, 202, 193, 100, 78, 208, 105, 87, 11, 53, 122, 16, 138, 180, 166, 194, 187, 4, 48, 134, 22, 190, 98, 95, 28, 186, 145, 109, 69, 18, 92, 130, 133, 52, 157, 81, 228, 200, 62, 238, 24, 68, 51, 102, 195, 61, 70, 50, 114, 213, 121, 26, 60, 234, 10, 65, 35, 42, 215, 128, 140, 58, 171, 217, 110, 126, 197, 77, 19, 156, 182, 175, 27, 17, 204, 173, 38, 147, 113, 1, 91, 74, 219, 101, 90, 150, 207, 123, 15, 57, 137, 106, 206, 88, 170, 45, 40, 93, 198, 47, 139, 199, 146, 118, 185, 41, 37, 46, 161, 76, 13, 216, 131, 83, 176, 230, 64, 7, 103, 205, 2, 152, 221, 25, 99, 136, 14, 233, 158, 141, 189, 119, 203, 3, 177, 89, 142, 148, 111, 226, 167, 67, 153, 31, 188, 143, 135, 104, 75, 192, 196, 168, 191, 97, 33, 54, 8, 229, 71, 144, 32, 79, 239, 112, 181, 227, 107, 240, 178, 115, 231, 172, 162, 237, 84, 155, 179, 72, 232, 164, 184, 160, 80, 159, 183]
best fitness: 5658
makespan: 180, total_wait: 2533, T_W_wait: 589
Gen: 324, best dna seq: [169, 21, 127, 212, 201, 43, 29, 210, 224, 116, 56, 9, 36, 96, 120, 209, 151, 223, 82, 94, 66, 154, 236, 225, 125, 23, 86, 129, 63, 132, 220, 108, 34, 59, 124, 174, 85, 55, 211, 12, 218, 222, 165, 235, 5, 149, 49, 39, 117, 73, 6, 30, 20, 163, 214, 44, 202, 193, 100, 78, 208, 105, 87, 11, 53, 122, 16, 138, 180, 166, 194, 187, 4, 48, 134, 22, 190, 98, 95, 28, 186, 145, 109, 69, 18, 92, 130, 133, 52, 157, 81, 228, 200, 62, 238, 24, 68, 51, 102, 195, 61, 70, 50, 114, 213, 121, 26, 60, 234, 10, 65, 35, 42, 215, 128, 140, 58, 171, 217, 110, 126, 197, 77, 19, 156, 182, 175, 27, 17, 204, 173, 38, 147, 113, 1, 91, 74, 219, 101, 90, 150, 207, 123, 15, 57, 137, 106, 206, 88, 170, 45, 40, 93, 198, 47, 139, 199, 146, 118, 185, 41, 37, 46, 161, 76, 13, 216, 131, 83, 176, 230, 64, 7, 103, 205, 2, 152, 221, 25, 99, 136, 14, 233, 158, 141, 189, 119, 203, 3, 177, 89, 142, 148, 111, 226, 167, 67, 153, 31, 188, 143, 135, 104, 75, 192, 196, 168, 191, 97, 33, 54, 8, 229, 71, 144, 32, 79, 239, 112, 181, 227, 107, 240, 178, 115, 231, 172, 162, 237, 84, 155, 179, 72, 232, 164, 184, 160, 80, 159, 183]
best fitness: 5658
makespan: 180, total_wait: 2533, T_W_wait: 589
Gen: 325, best dna seq: [169, 21, 127, 212, 201, 43, 29, 210, 224, 116, 56, 9, 36, 96, 120, 209, 151, 223, 82, 94, 66, 154, 236, 225, 125, 23, 86, 129, 63, 132, 220, 108, 34, 59, 124, 174, 85, 55, 211, 12, 218, 222, 165, 235, 5, 149, 49, 39, 117, 73, 6, 30, 20, 163, 214, 44, 202, 193, 100, 78, 208, 105, 87, 11, 53, 122, 16, 138, 180, 166, 194, 187, 4, 48, 134, 22, 190, 98, 95, 28, 186, 145, 109, 69, 18, 92, 130, 133, 52, 157, 81, 228, 200, 62, 238, 24, 68, 51, 102, 195, 61, 70, 50, 114, 213, 121, 26, 60, 234, 10, 65, 35, 42, 215, 128, 140, 58, 171, 217, 110, 126, 197, 77, 19, 156, 182, 175, 27, 17, 204, 173, 38, 147, 113, 1, 91, 74, 219, 101, 90, 150, 207, 123, 15, 57, 137, 106, 206, 88, 170, 45, 40, 93, 198, 47, 139, 199, 146, 118, 185, 41, 37, 46, 161, 76, 13, 216, 131, 83, 176, 230, 64, 7, 103, 205, 2, 152, 221, 25, 99, 136, 14, 233, 158, 141, 189, 119, 203, 3, 177, 89, 142, 148, 111, 226, 167, 67, 153, 31, 188, 143, 135, 104, 75, 192, 196, 168, 191, 97, 33, 54, 8, 229, 71, 144, 32, 79, 239, 112, 181, 227, 107, 240, 178, 115, 231, 172, 162, 237, 84, 155, 179, 72, 232, 164, 184, 160, 80, 159, 183]
best fitness: 5658
makespan: 180, total_wait: 2533, T_W_wait: 589
```

## 实验

遗传算法算例计算10次，取其中最优解

种群规模： 50, 100, **<span style="color:red;">150</span>**, 200;

交叉概率：0.5, **<span style="color:red;">0.8</span>**, 1.0

变异概率:0.7, 0.8, 0.9, **<span style="color:red;">1.0</span>**

等待时间阈值设置为15分钟

单项目等待时间超过阈值惩罚成本参数为5





