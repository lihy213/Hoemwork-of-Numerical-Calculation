#!/usr/bin/python3
# -*- coding : utf-8 -*-
# @Project : Pycharm Project
# @Subproject  : conjugte_gradient_exercise
# @Author : howyoung
# @Time :  2022/10/26 11:15

import numpy as np
import matplotlib.pyplot as plt


class ConjugateGradient(object):
    """
    共轭梯度法
    A为系数矩阵（在这里为对称正定矩阵）
    b为方程的右端项
    n为系数矩阵的维数
    precision为设置的精度
    x0为给定的迭代初始向量
    iter为迭代次数
    """

    def __init__(self, A, b, n, precision, x0):
        self.A = np.array(A)
        self.b = b
        self.n = n
        self.precision = precision
        self.x0 = np.array(x0)
        self.x = []
        self.x.append(list(x0))   # 需要首先让x中有初始值
        self.iter = 0
        self.f = open('共轭梯度求解结果.log', 'a', encoding='utf-8')

    def norm_2(self, vector):
        """
        求解2范数
        :param vector: vector refers one specific vector(here refers to
        """
        sum_norm_2 = [sum(i ** 2 for i in vector)]

        return np.sqrt(sum_norm_2)

    def cal(self):
        """
        按照共轭梯度法流程求解
        """
        r0 = self.b - np.dot(self.A, self.x0)   # 为n*1的一维矩阵
        d0 = r0
        d_acr = d0    # 初始d，r准确值赋值
        r_acr = r0
        result_residual = []
        while True:
            # 计算d_k_T与A的乘积
            temp = np.dot(d_acr.T, self.A)
            # 计算alpha_k
            alpha = np.dot(r_acr.T, r_acr) / np.dot(temp, d_acr)
            # 计算x_k+1
            x = np.array(self.x[-1]) + np.dot(alpha, d_acr)
            # 计算残差r_k+1
            r = self.b - np.dot(self.A, x)
            # 将计算出的x加入到向量x中
            self.x.append(list(x))
            result_residual.append(np.log10(self.norm_2(r)))
            self.iter += 1
            # 判断残差与给定值的大小
            if (self.norm_2(r) <= self.precision) or (self.iter + 1 == self.n):
                break
            else:
                # 若不满足，则继续求得beta和d_k+1
                beta = (self.norm_2(r) ** 2) / (self.norm_2(r_acr) ** 2)
                d = r + beta * d_acr
            r_acr = r
            d_acr = d

        return np.array(result_residual), self.iter

    def picture_show_error(self):
        """绘制error图像"""

        # 调用cal()函数
        a = []
        for i in self.cal()[0]:
            a.append(i[0])
        hor_line = np.linspace(1, self.iter, len(self.x) - 1)
        print("="*20, file=self.f)
        print("\n第", self.iter, "计算的近似结果为:", self.x[-1], file=self.f)
        print("\n迭代次数:", self.iter, file=self.f)

        plt.plot(hor_line, np.array(a), marker='*', c='r', linestyle='--', label='error-iterations')
        plt.title("error", fontsize=16, c='r')
        plt.xlabel("iterations", fontsize=16)
        plt.ylabel("$\mathrm{log}$$_\mathrm{10}$$\mathrm{(error)}$", fontsize=16)
        plt.tick_params(axis='both', labelsize=14)
        plt.grid(True)
        plt.savefig('n=' + str(self.n) + '.jpg', bbox_inches='tight', dpi=640)
        plt.show()
        plt.close()


if __name__ == '__main__':
    for n in [100, 200, 400]:
        # 构造矩阵A
        A = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                if i == 0:
                    A[i][0] = -2
                    A[i][1] = 1
                elif i == n - 1:
                    A[i][n - 1] = -2
                    A[i][n - 2] = 1
                elif i == j:
                    A[i][j] = -2
                    A[i][j - 1] = 1
                    A[i][j + 1] = 1
        # 构造向量b
        b = [0] * n
        b[-1] = -1
        b[0] = -1
        # 初始向量赋值[0 0 0 0 0 ……]
        x0 = [0] * n
        # 计算精度
        precision = 1e-10

        start = ConjugateGradient(A, b, n, precision, x0)
        start.picture_show_error()
