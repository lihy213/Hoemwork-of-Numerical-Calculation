#!/usr/bin/python3
# -*- coding : utf-8 -*-
# @Project : Iteration
# @Subproject  : erfenfa
# @Author : howyoung
# @Software : PyCharm
# @Time :  2022/11/16 18:26

import matplotlib.pyplot as plt


class ErFenFa(object):

    def __init__(self, a, b):
        """
        定义初始变量
        :param a: 区间左端点
        :param b: 区间右端点
        """
        self.a0 = a
        self.b0 = b
        self.x0 = 0.5 * (self.a0 + self.b0)

        # 定义误差精度，epsilon1、epsilon2分别为方程值和区间长度对应的误差精度
        self.epsilon1 = 1e-8
        self.epsilon2 = 1e-8
        self.k = 0
        self.list_fx = []
        self.list_x = []
        self.list_delta = []

    def func(self, x):
        """
        求解方程
        :param x:自变量，此处由cal_fx传入
        :return: 方程的值f(x)
        """
        fx = pow(x, 6) - 5 * pow(x, 5) + 3 * pow(x, 4) + pow(x, 3) - 7 * pow(x, 2) + 7 * x - 20
        return fx

    def cal_fx(self):
        """
        函数主体部分，用于判断收敛性，并变换区间
        :return: 满足收敛条件的方程的近似根及一些残差
        """
        x_k = self.x0
        a_k = self.a0
        b_k = self.b0
        f = open('result.log', 'w', encoding='utf-8')
        print('=' * 10, '迭代结果', '=' * 10, file=f)

        while self.k <= 100:
            self.k += 1
            fx_k = self.func(x_k)
            delta_x = abs(b_k - a_k)
            self.list_fx.append(fx_k)
            self.list_delta.append(delta_x)
            self.list_x.append(x_k)
            # 输出结果到文本文件
            print('\n第', self.k, '次迭代', file=f)
            print('fx误差：', fx_k, end='   ', file=f)
            print('区间长度残差：', delta_x, end='   ', file=f)
            print('方程的近似根：', x_k, file=f)

            if abs(fx_k) < self.epsilon1 or delta_x < self.epsilon2:
                print('计算收敛！', file=f)
                return x_k, self.list_x, self.list_fx, self.list_delta
            else:
                if self.func(a_k)*self.func(x_k) < 0:
                    a_kp1 = a_k
                    b_kp1 = x_k
                    x_k = 0.5*(a_kp1 + b_kp1)
                    # 此处比较关键，需要获取b的最新值
                    b_k = b_kp1
                else:
                    a_kp1 = x_k
                    b_kp1 = b_k
                    x_k = 0.5 * (a_kp1 + b_kp1)
                    # 此处比较关键，需要获取a的最新值
                    a_k = a_kp1

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
        plt.annotate(text=(self.k, vert[-1]), xy=(self.k, vert[-1]), xytext=(0.5*self.k, max(vert)),
                     weight='bold', color='r', arrowprops=dict(arrowstyle='->', connectionstyle='arc3', color='c'),
                     bbox=dict(boxstyle='round,pad=0.5', fc='yellow', ec='k', lw=1, alpha=0.4))
        plt.scatter(self.k, vert[-1], marker='v', s=50, c='r')
        plt.grid(True)
        plt.savefig('Residual_fx.png', bbox_inches='tight')
        plt.show()

    def plot_delta(self):
        """
        绘制关于区间长度的误差图
        """
        horizon = [hor for hor in range(1, self.k+1)]
        vert = self.list_delta

        plt.title('Residual_delta', fontsize=20, c='r')
        plt.xlabel('Iterations', fontsize=16)
        plt.ylabel('Abs(b_k-a_k)', fontsize=16)
        plt.plot(horizon, vert, linewidth=2)
        plt.annotate(text=(self.k, vert[-1]), xy=(self.k, vert[-1]), xytext=(0.5*self.k, max(vert)),
                     weight='bold', color='r', arrowprops=dict(arrowstyle='->', connectionstyle='arc3', color='c'),
                     bbox=dict(boxstyle='round,pad=0.5', fc='yellow', ec='k', lw=1, alpha=0.4))
        plt.scatter(self.k, vert[-1], marker='h', s=50, c='r')
        plt.grid(True)
        plt.savefig('Residual_delta.png', bbox_inches='tight')
        plt.show()


if __name__ == '__main__':
    cal = ErFenFa(-1, 5)
    result = cal.cal_fx()
    cal.plot_fx()
    cal.plot_delta()
