#!/usr/bin/python3
# -*- coding : utf-8 -*-
# @Project : project2
# @Subproject  : SquareLeastMethod
# @Author : howyoung
# @Software : PyCharm
# @Time :  2022/10/31 8:49


import numpy as np
import GS as gs
import Condugate_gradient_method as cg


class SquareLeastMethod(object):

    def __init__(self, X, Y, n):
        """
        给定初始条件
        :param X: 题目中的X
        :param Y: 题目中的Y
        :param n:求解的最小二乘多项式的最高次次数
        """
        self.X = X
        self.Y = Y
        self.n = n + 1
        self.A = np.ones((self.n, self.n))
        self.b = np.ones(self.n)
        self.f = open('求解结果.log', 'w', encoding='utf-8')

    def cal__coeff_matrix(self):
        """
        求解系数矩阵和方程右端项
        """
        for i in range(self.n):
            for j in range(i, self.n):
                # 求解对角线元素
                if i == j:
                    self.A[i, j] = sum(pow(g, 2*i) for g in self.X)
                else:
                    self.A[i, j] = sum(pow(m, j+i) for m in self.X)
                    self.A[j, i] = self.A[i, j]
            # 求解右端项
            X_new = [pow(k, i) for k in self.X]
            self.b[i] = np.dot(self.Y, X_new)

        print('系数矩阵为:\n', self.A, file=self.f)
        print('\n方程右端项为:\n', self.b, file=self.f)

    def switch(self, switch):
        """
        调用数值迭代方法
        :param switch: 迭代方法，type: str; gs表示Gauss-Sieder方法；cg表示共轭梯度法
        """
        if switch == 'cg':
            # 调用共轭梯度法程序
            x0 = np.zeros(5)
            start = cg.ConjugateGradient(self.A, self.b, self.n, 1e-5, x0)
            start.picture_show_error()
        elif switch == 'gs':
            # 调用Gauss-Sieder程序
            X = gs.GaussMethod(self.A, self.b, 100000, 1e-5)
            X.plot_residuals()
            var = X.Jacobi()
            print("==========Gauss-Sieder计算结果===========", file=self.f)
            print('\n求解结果为:\n', var[1], file=self.f)
            print('\n',self.n-1, '次多项式为:\n', 'Pn(x)=',var[1][0],'+',
                  var[1][1],'*x','+',var[1][2],'*x^2','+',
                  var[1][3],'*x^3','+',var[1][4],'*x^4', file=self.f)


if __name__ == '__main__':
    # 书5.1计算实习例题
    X = np.array(np.linspace(0.1, 0.9, 9))
    Y = np.array([5.1234, 5.3057, 5.5687,
                  5.9375, 6.4370, 7.0978,
                  7.9493, 9.0253, 10.3627])
    SLM = SquareLeastMethod(X, Y, n=4)
    SLM.cal__coeff_matrix()
    # 通过选择switch从而选择对应的迭代方法
    # gs表示Gauss-Sieder方法；cg表示共轭梯度法
    SLM.switch('cg')
