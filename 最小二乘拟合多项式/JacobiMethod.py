#!/usr/bin/python3
# -*- coding : utf-8 -*-
# @Project : project2
# @Subproject  : JacobiMethod
# @Author : howyoung
# @Software : PyCharm
# @Time :  2022/11/1 19:04

import numpy as np


def Jacobi(A, b, k, tol):
    """
    雅可比迭代
    :param A:方程系数矩阵
    :param b: 方程右端项
    :param k: 给定的最大迭代次数
    :param tol: 给定计算残差
    :return: 方程组的解
    """
    f = open('Jacobi.log', 'a', encoding='utf-8')
    n = A.shape[1]
    # 输出对角为1的矩阵D
    D = np.eye(n)
    # 将D的对角线元素1置换为A的对角线元素
    D[np.arange(n), np.arange(n)] = A[np.arange(n), np.arange(n)]
    LU = A - D
    X = np.zeros(n)

    # 迭代k次
    for i in range(k):
        # 将上一次求解的X赋值给最新的X1
        X1 = X.copy()
        # 求矩阵D的逆矩阵
        D_inv = np.linalg.inv(D)
        # 按照数值求解方法解出X
        X = -np.dot(np.dot(D_inv, LU), X) + np.dot(D_inv, b)

        print('\n第', i, '次迭代：', X, file=f)
        print('残差值为：', np.abs(X - X1), file=f)
        print('当前最大残差为：', np.max(np.abs(X - X1)), file=f)

        if np.max(np.abs(X - X1)) < tol:
            print('计算收敛!', file=f)
            break

    return X


A = np.array([
    [3, -1, 0, 0],
    [-2, 6, -1, 0],
    [0, -2, 6, -1],
    [0, 0, -2, 7]
])
b = np.array([3, 4, 5, -3])
X = Jacobi(A, b, 10000, 1e-5)
