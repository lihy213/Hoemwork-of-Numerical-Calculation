#!/usr/bin/python3
# -*- coding : utf-8 -*-
# @Project : Iteration
# @Subproject  : easy_iteration
# @Author : howyoung
# @Software : PyCharm
# @Time :  2022/11/16 21:38

import matplotlib.pyplot as plt


class EasyIteration(object):

    def __init__(self, a, b):
        self.a0 = a
        self.b0 = b
        self.x0 = 0.5 * (self.a0 + self.b0)
        self.k = 0
        self.epsilon1 = 1e-8
        self.epsilon2 = 1e-8
        self.list_fx = []
        self.list_x = []
        self.list_delta = []

    def func(self, x):
        """
        求解方程
        :param x:自变量，此处由cal_fx传入
        :return: 方程的值f(x), 迭代phi(x)
        """
        fx = pow(x, 6) - 5 * pow(x, 5) + 3 * pow(x, 4) + pow(x, 3) - 7 * pow(x, 2) + 7 * x - 20
        phi_x = pow(5 * pow(x, 5) - 3 * pow(x, 4) - pow(x, 3) + 7 * pow(x, 2) - 7 * x + 20, 1/6)
        return fx, phi_x

    def cal(self):
        """
        简单迭代法收敛性判断
        :return: 满足收敛条件的方程的近似根及一些残差
        """
        # 首先使用一步二分法求出一个初始值x1；此处未进行收敛性判定，因为假定第一次不可能收敛，也符合实际情况
        if self.func(self.a0)[0] * self.func(self.x0)[0] < 0:
            a_kp1 = self.a0
            b_kp1 = self.x0
            x_k = 0.5 * (a_kp1 + b_kp1)
        else:
            a_kp1 = self.x0
            b_kp1 = self.b0
            x_k = 0.5 * (a_kp1 + b_kp1)
        # x_k = self.x0
        # 得到初始delta_x
        delta_x = (self.x0 + x_k)/2
        # delta_x = self.epsilon2
        f = open('result_easy.log', 'w', encoding='utf-8')
        print('=' * 10, '迭代结果', '=' * 10, file=f)

        while True:
            self.k += 1
            fx_k = self.func(x_k)[0]
            self.list_fx.append(fx_k)
            self.list_x.append(x_k)
            self.list_delta.append(delta_x)

            # 输出结果到文本文件
            print('\n第', self.k, '次迭代', file=f)
            print('fx误差：', fx_k, end='   ', file=f)
            print('区间长度残差：', delta_x, end='   ', file=f)
            print('方程的近似根：', x_k, file=f)

            if abs(fx_k) < self.epsilon1 or delta_x < self.epsilon2:
                print('计算收敛！', file=f)
                return x_k, self.list_x, self.list_fx, self.list_delta
            else:
                x_next = self.func(x_k)[1]
                delta_x = abs(x_next - x_k)
                x_k = x_next

    def plot_fx(self):
        """
        绘制关于方程值的误差图
        """
        horizon = [hor for hor in range(1, self.k+1)]
        vert = self.list_fx

        plt.title('Residual_fx', fontsize=20, c='r')
        plt.xlabel('Iterations', fontsize=16)
        plt.ylabel('fx', fontsize=16)
        plt.plot(horizon, vert, linewidth=2)
        plt.annotate(text=(self.k, vert[-1]), xy=(self.k, vert[-1]), xytext=(0.5*self.k, min(vert)),
                     weight='bold', color='r', arrowprops=dict(arrowstyle='->', connectionstyle='arc3', color='c'),
                     bbox=dict(boxstyle='round,pad=0.5', fc='yellow', ec='k', lw=1, alpha=0.4))
        plt.scatter(self.k, vert[-1], marker='v', s=50, c='r')
        plt.grid(True)
        plt.savefig('Residual_fx_easy.png', bbox_inches='tight')
        plt.show()

    def plot_delta(self):
        """
        绘制相邻两个近似根的误差图
        """
        horizon = [hor for hor in range(1, self.k+1)]
        vert = self.list_delta

        plt.title('Residual_delta', fontsize=20, c='r')
        plt.xlabel('Iterations', fontsize=16)
        plt.ylabel('Abs(x_k+1-x_k)', fontsize=16)
        plt.plot(horizon, vert, linewidth=2)
        plt.annotate(text=(self.k, vert[-1]), xy=(self.k, vert[-1]), xytext=(0.5*self.k, max(vert)),
                     weight='bold', color='r', arrowprops=dict(arrowstyle='->', connectionstyle='arc3', color='c'),
                     bbox=dict(boxstyle='round,pad=0.5', fc='yellow', ec='k', lw=1, alpha=0.4))
        plt.scatter(self.k, vert[-1], marker='h', s=50, c='r')
        plt.grid(True)
        plt.savefig('Residual_delta_easy.png', bbox_inches='tight')
        plt.show()


if __name__ == '__main__':
    cal = EasyIteration(-1, 5)
    result = cal.cal()
    cal.plot_fx()
    cal.plot_delta()
