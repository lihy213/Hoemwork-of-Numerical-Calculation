#!/usr/bin/python3
# -*- coding : utf-8 -*-
# @Project : Iteration
# @Subproject  : NewtonMethod
# @Author : howyoung
# @Software : PyCharm
# @Time :  2022/11/17 10:42

import sympy as sp
import numpy as np
import matplotlib.pyplot as plt


class NewTon(object):

    def __init__(self, a, b):
        """
        一些初始条件
        :param a: 区间左端点
        :param b: 区间右端点
        """
        self.a0 = a
        self.b0 = b
        self.x0 = 0.5 * (self.a0 + self.b0)
        self.k = 0
        self.epsilon1 = 1e-8
        self.epsilon2 = 1e-8
        self.list_fx = []
        self.list_x = []
        self.list_delta = []
        self.f = open('result_Newton.log', 'w', encoding='utf-8')

    def func(self):
        """
        求解方程
        :return: 方程的值f(x), 方程的一阶导数D_fx和二阶导数D2_fx
        """
        x = sp.symbols('x')
        fx = pow(x, 6) - 5 * pow(x, 5) + 3 * pow(x, 4) + pow(x, 3) - 7 * pow(x, 2) + 7 * x - 20
        D_fx = sp.diff(fx, x)
        D2_fx = sp.diff(fx, x, 2)
        return fx, D_fx, D2_fx

    def cal_fx(self, x_k):
        # 求解fx
        fx_value = self.func()[0].evalf(subs={'x': x_k})

        return fx_value

    def cal_Dfx(self, x_k):
        # 求解D_fx
        Dfx_value = self.func()[1].evalf(subs={'x': x_k})
        return Dfx_value

    def main_cal(self):
        """
        主函数，用于判断是否达到拟定收敛精度，并输出相应结果
        :return: 方程的近似根，方程的值，相邻两步迭代近似根的误差
        """
        x_k = self.x0
        delta_x = self.epsilon2

        print('=' * 10, '迭代结果', '=' * 10, file=self.f)

        while True:
            self.k += 1
            self.list_fx.append(self.cal_fx(x_k))
            x_next = x_k - self.cal_fx(x_k)/self.cal_Dfx(x_k)
            self.list_x.append(x_k)
            self.list_delta.append(delta_x)

            # 输出结果到文本文件
            print('\n第', self.k, '次迭代', file=self.f)
            print('fx误差：', self.cal_fx(x_k), end='   ', file=self.f)
            print('区间长度残差：', delta_x, end='   ', file=self.f)
            print('方程的近似根：', x_k, file=self.f)

            if delta_x < self.epsilon2 or abs(self.cal_fx(x_k)) < self.epsilon1:
                print('计算收敛！\n', file=self.f)
                return x_k, self.list_x, self.list_fx, self.list_delta
            else:
                delta_x = abs(x_next - x_k)
                x_k = x_next
    """绘图"""
    def cal_D2fx(self):
        """
        计算D2_fx，并找出其中相邻两个值异号的点
        :return: D2_fx
        """
        list_D2fx = [0]
        for i in np.linspace(-1, 5, 1000):
            # a = self.func(i)
            a = self.func()[2].evalf(subs={'x': i})
            # list_D2fx.append(a[2])
            list_D2fx.append(a)
            if list_D2fx[-1]*list_D2fx[-2] < 0:
                # 寻找并标记异号的D2_fx
                plt.scatter(i, a, marker='v', s=50, c='r')
                plt.annotate(text=(i, a), xy=(i, a), xytext=(i, abs(1000*list_D2fx[-1])),
                             weight='bold', color='r',
                             arrowprops=dict(arrowstyle='->', connectionstyle='arc3', color='c'),
                             bbox=dict(boxstyle='round,pad=0.5', fc='yellow', ec='k', lw=1, alpha=0.4))

        return list_D2fx

    def plot_D2fx(self):
        """
        绘制fx二阶导数D2_fx随x的变化关系图
        """
        vert = self.cal_D2fx()
        # 判断全局收敛性，是否满足方程二阶导数不变号的性质
        if max(vert) * min(vert) < 0:
            print('不满足二阶导数不变号的条件...<<<全局不收敛...', file=self.f)
        # 绘图主体部分
        plt.title('Variation of D2_fx with x', fontsize=20, c='r', fontweight='bold')
        plt.xlabel('Value of x', fontsize=16)
        plt.ylabel('Value of D2_fx', fontsize=16)
        plt.text(-1, max(vert)/1.2, s='D2_fx表示函数fx的二阶导数', fontproperties='SimHei',
                 bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 5})
        plt.plot(np.linspace(-1, 5, 1000), np.array(vert[1:]))
        plt.grid(True)
        plt.savefig('D2fx_Newton.png', bbox_inches='tight')
        plt.show()

    def plot_fx(self):
        """
        绘制关于方程值的误差图
        """
        horizon = [hor for hor in range(1, self.k+1)]
        vert = self.list_fx

        plt.title('Residual_fx', fontsize=20, c='r', fontweight='bold')
        plt.xlabel('Iterations', fontsize=16)
        plt.ylabel('fx', fontsize=16)
        plt.plot(horizon, vert, linewidth=2)
        plt.annotate(text=(self.k, vert[-1]), xy=(self.k, vert[-1]), xytext=(0.5*self.k, max(vert)),
                     weight='bold', color='r', arrowprops=dict(arrowstyle='->', connectionstyle='arc3', color='c'),
                     bbox=dict(boxstyle='round,pad=0.5', fc='yellow', ec='k', lw=1, alpha=0.4))
        plt.scatter(self.k, vert[-1], marker='v', s=50, c='r')
        plt.grid(True)
        plt.savefig('Residual_fx_Newton.png', bbox_inches='tight')
        plt.show()

    def plot_delta(self):
        """
        绘制相邻两个近似根的误差图
        """
        horizon = [hor for hor in range(1, self.k+1)]
        vert = self.list_delta

        plt.title('Residual_delta', fontsize=20, c='r', fontweight='bold')
        plt.xlabel('Iterations', fontsize=16)
        plt.ylabel('Abs(x_k+1-x_k)', fontsize=16)
        plt.plot(horizon, vert, linewidth=2)
        plt.annotate(text=(self.k, vert[-1]), xy=(self.k, vert[-1]), xytext=(0.5*self.k, max(vert)),
                     weight='bold', color='r', arrowprops=dict(arrowstyle='->', connectionstyle='arc3', color='c'),
                     bbox=dict(boxstyle='round,pad=0.5', fc='yellow', ec='k', lw=1, alpha=0.4))
        plt.scatter(self.k, vert[-1], marker='h', s=50, c='r')
        plt.grid(True)
        plt.savefig('Residual_delta_Newton.png', bbox_inches='tight')
        plt.show()

    def convergency(self):
        """
        判断全局收敛性
        """
        for i in np.linspace(-1, 5, 1000):
            var_Dfx = self.cal_Dfx(i)
            if var_Dfx == 0:
                print('满足一阶导数不为零的条件...<<<全局不收敛...', file=self.f)
        else:
            print('满足一阶导数不为零的条件...', file=self.f)
        var_fx1 = self.cal_fx(self.a0)
        var_fx2 = self.cal_fx(self.b0)
        if var_fx1 * var_fx2 >= 0:
            print('不满足区间端点方程值乘积小于零的条件...<<<全局不收敛...', file=self.f)
        else:
            print('满足区间端点方程值乘积小于零的条件...', file=self.f)


if __name__ == '__main__':
    cal = NewTon(-1, 5)
    cal.main_cal()
    cal.cal_D2fx()
    cal.plot_D2fx()
    cal.convergency()
    cal.plot_fx()
    cal.plot_delta()
