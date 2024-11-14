#!/usr/bin/python3
# -*- coding : utf-8 -*-
# @Project : project2
# @Subproject  : Gauss-SiederMethod
# @Author : howyoung
# @Software : PyCharm
# @Time :  2022/11/1 19:05

import numpy as np
import matplotlib.pylab as plt


class GaussMethod(object):

    def __init__(self, A, b, iter, tol):
        """定义初始值"""
        self.A = A
        self.b = b
        self.iter = iter
        self.tol = tol
        self.residuals = []

    def Jacobi(self):
        """Jacobi迭代求解"""
        f = open("Gauss-Sieder_residual.log", 'w', encoding='utf-8')

        n = self.A.shape[1]
        D = np.eye(n)
        D[np.arange(n), np.arange(n)] = self.A[np.arange(n), np.arange(n)]
        LU = self.A - D
        # 得到LU的上、下三角矩阵
        L = np.tril(LU)
        U = np.triu(LU)
        D_L = D + L
        X = np.zeros(n)

        # 最多迭代iter次
        for i in range(self.iter):
            X1 = X.copy()
            D_L_inv = np.linalg.inv(D_L)
            X = -np.dot(np.dot(D_L_inv, U), X) + np.dot(D_L_inv, self.b)

            print('\n第', i, '次迭代：', X, file=f)
            print('残差值为：', np.abs(X - X1), file=f)
            print('当前最大残差为：', np.max(np.abs(X - X1)), file=f)

            self.residuals.append(np.max(np.abs(X - X1)))

            if np.max(np.abs(X - X1)) < self.tol:
                print('计算收敛!', file=f)
                break

        f.close()

        return i, X

    def plot_residuals(self):
        """绘制残差曲线"""
        var = self.Jacobi()[0]
        horizon = np.linspace(0, var, var + 1)
        verticle = np.log10(self.residuals)

        plt.title("Residual", fontsize=20)
        plt.plot(horizon, verticle, marker='*', linestyle=':', c='r', label='Residual', linewidth=2)
        plt.xlabel("Iterations", fontsize=16)
        plt.ylabel("log$_\mathit{10}$(residual)", fontsize=16)
        plt.grid(True)
        plt.tick_params(axis='both')
        # 显示收敛点或迭代结束点
        plt.scatter(var, verticle[-1], s=50, c='b')
        plt.text(var, verticle[-1], (var, verticle[-1]), c='b')
        plt.savefig('Gauss-Sieder.png', bbox_inches='tight')
        plt.show()
