# -*- coding: UTF-8 -*-
import numpy as np
import math
import os

path = 'data/1.2'
dirs = os.listdir(path)
data = []
print(dirs)
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

data_2 = data[7:9]
f_0_1 = ((np.mean(data_2[0][:, 3], 0)) / 1000 +
         (np.mean(data_2[1][:, 3], 0)) / 1000) / 2
print(f"f_0_1 = {f_0_1}")

# ----------------------------1.2.3----------------------------

f_0_2_test = ((np.mean(data[10][:, 3], 0) / 1000) +
              (np.mean(data[12][:, 3], 0) / 1000)) / 2
f_0_3_test = ((np.mean(data[9][:, 3], 0) / 1000) +
              (np.mean(data[11][:, 3], 0) / 1000)) / 2
B = [f_0_1, f_0_2_test, f_0_3_test]
B_r = np.std(B, ddof=1)
print(f"B_r = {B_r}")
