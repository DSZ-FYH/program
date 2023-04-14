# *-* coding: utf-8 *-*
import numpy as np
import matplotlib.pyplot as plt
import math

deg = 180/math.pi

data_save = np.loadtxt('data_save.txt')
GPS = np.loadtxt('GPS.txt')
sins_save = np.loadtxt('sins_save.txt')
P_save = np.loadtxt('P_save.txt')

# print(data_save[188945])

fig = plt.figure('3D滤波结果')
ax = fig.add_subplot(111, projection='3d')
ax.plot(data_save[:,0]*deg,data_save[:,1]*deg,data_save[:,2]*deg)
ax.plot(GPS[:,2]*deg,GPS[:,3]*deg,GPS[:,4]*deg)
ax.plot(sins_save[:,0]*deg,sins_save[:,1]*deg,sins_save[:,2]*deg)
ax.set_xlabel('纬度/°')
ax.set_ylabel('经度/°')
ax.set_zlabel('高度/h')
ax.legend(['滤波轨迹','GPS原始数据','纯惯性轨迹'])
ax.grid(True)
plt.show()

fig = plt.figure('2D滤波结果')
plt.plot(data_save[:,0]*deg,data_save[:,1]*deg)
plt.plot(GPS[:,2]*deg,GPS[:,3]*deg)
plt.plot(sins_save[:,0]*deg,sins_save[:,1]*deg)
plt.title('2D滤波结果')
plt.xlabel('纬度/°')
plt.ylabel('经度/°')
plt.legend(['滤波轨迹','GPS原始数据','纯惯性轨迹'])
plt.grid(True)
plt.show()

fig = plt.figure('位置滤波结果')
plt.subplot(3,1,1)
plt.plot(data_save[:,0]*deg)
plt.plot(GPS[:,2]*deg)
plt.title('纬度滤波结果')
plt.xlabel('时间')
plt.ylabel('纬度/°')
plt.legend(['滤波轨迹','GPS原始数据'])
plt.grid(True)

