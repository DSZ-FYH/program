# -*- coding: UTF-8 -*-
import numpy as np
import math
import os

path = 'data/1.2'
dirs = os.listdir(path)
data = []
for file in dirs:
    data.append(np.loadtxt(os.path.join(path, file)))

# ----------------------------1.2.1----------------------------
data_1 = data[0:7]

g = 9.7803267714
X = []
Y = []
for i in range(0, len(data_1)):
    X.append(g * math.cos(math.pi * i / 3))
    Y.append(np.mean(data_1[i][:, 3], 0) / 1000)

fit = np.polyfit(X, Y, deg=1)

K_A = fit[0]
F_A0 = fit[1] / 7
print(f"K_A = {K_A}")
print(f"F_A0 = {F_A0}")

# ----------------------------1.2.2----------------------------
data_a_1 = np.loadtxt("data/1.2/1.2-2 0度.txt")
data_b_1 = np.loadtxt("data/1.2/1.2-2 180度.txt")

f_0_1 = ((np.mean(data_a_1[:, 3], 0) / 1000) +
         (np.mean(data_b_1[:, 3], 0) / 1000)) / 2
print(f"f_0 = {f_0_1}")

# ----------------------------1.2.3----------------------------
data_a_2 = np.loadtxt("data/1.2/1.2-3 0度 二组.txt")
data_a_3 = np.loadtxt("data/1.2/1.2-3 0度 三组.txt")
data_b_2 = np.loadtxt("data/1.2/1.2-3 180度 二组.txt")
data_b_3 = np.loadtxt("data/1.2/1.2-3 180度 三组.txt")
f_0_2 = ((np.mean(data_a_2[:, 3], 0) / 1000) +
         (np.mean(data_b_2[:, 3], 0) / 1000)) / 2
f_0_3 = ((np.mean(data_a_3[:, 3], 0) / 1000) +
         (np.mean(data_b_3[:, 3], 0) / 1000)) / 2
B = [f_0_1, f_0_2, f_0_3]
B_r = np.std(B, ddof=1)
print(f"B_r = {B_r}")
