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

fig = plt.figure('3D滤波结果',figsize=(8,8))
ax = fig.add_subplot(111, projection='3d')
ax.plot(data_save[:,0]*deg,data_save[:,1]*deg,data_save[:,2]*deg)
ax.plot(GPS[:,2]*deg,GPS[:,3]*deg,GPS[:,4]*deg)
ax.plot(sins_save[:,0]*deg,sins_save[:,1]*deg,sins_save[:,2]*deg)
plt.rcParams['font.sans-serif']=['SimHei']
plt.rcParams['axes.unicode_minus']=False
ax.set_xlabel('纬度/°')
ax.set_ylabel('经度/°')
ax.set_zlabel('高度/h')
ax.legend(['滤波轨迹','GPS原始数据','纯惯性轨迹'])
ax.grid(True)
plt.show()

fig = plt.figure('2D滤波结果',figsize=(8,8))
plt.plot(data_save[:,0]*deg,data_save[:,1]*deg)
plt.plot(GPS[:,2]*deg,GPS[:,3]*deg)
plt.plot(sins_save[:,0]*deg,sins_save[:,1]*deg)
plt.title('2D滤波结果')
plt.rcParams['font.sans-serif']=['SimHei']
plt.rcParams['axes.unicode_minus']=False
plt.xlabel('纬度/°')
plt.ylabel('经度/°')
plt.legend(['滤波轨迹','GPS原始数据','纯惯性轨迹'])
plt.grid(True)
plt.show()

fig = plt.figure('位置滤波结果',figsize=(16,9))
plt.subplot(311)
plt.plot(data_save[:,0]*deg)
plt.plot(GPS[:,2]*deg)
plt.title('纬度滤波结果')
plt.rcParams['font.sans-serif']=['SimHei']
plt.rcParams['axes.unicode_minus']=False
plt.xlabel('时间')
plt.ylabel('纬度/°')
plt.legend(['滤波轨迹','GPS原始数据'])
plt.grid(True)

plt.subplot(312)
plt.plot(data_save[:,1]*deg)
plt.plot(GPS[:,3]*deg)
plt.title('经度滤波结果')
plt.rcParams['font.sans-serif']=['SimHei']
plt.rcParams['axes.unicode_minus']=False
plt.xlabel('时间')
plt.ylabel('经度/°')
plt.legend(['滤波轨迹','GPS原始数据'])
plt.grid(True)

plt.subplot(313)
plt.plot(data_save[:,2]*deg)
plt.plot(GPS[:,4]*deg)
plt.title('高度滤波结果')
plt.rcParams['font.sans-serif']=['SimHei']
plt.rcParams['axes.unicode_minus']=False
plt.xlabel('时间')
plt.ylabel('高度/h')
plt.legend(['滤波轨迹','GPS原始数据'])
plt.grid(True)
plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.2, hspace=0.4)
plt.show()

fig = plt.figure('姿态滤波结果',figsize=(16,9))
plt.subplot(311)
plt.plot(data_save[:,3]*deg)
plt.title('俯仰滤波结果')
plt.rcParams['font.sans-serif']=['SimHei']
plt.rcParams['axes.unicode_minus']=False
# plt.xlabel('时间')
plt.ylabel('俯仰/°')
plt.grid(True)

plt.subplot(312)
plt.plot(data_save[:,4]*deg)
plt.title('横滚滤波结果')
plt.rcParams['font.sans-serif']=['SimHei']
plt.rcParams['axes.unicode_minus']=False
# plt.xlabel('时间')
plt.ylabel('横滚/°')
plt.grid(True)

plt.subplot(313)
plt.plot(data_save[:,5]*deg)
plt.title('偏航滤波结果')
plt.rcParams['font.sans-serif']=['SimHei']
plt.rcParams['axes.unicode_minus']=False
# plt.xlabel('时间')
plt.ylabel('偏航/°')
plt.grid(True)
plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.2, hspace=0.4)
plt.show()

fig = plt.figure('速度滤波结果',figsize=(16,9))
plt.subplot(311)
plt.plot(data_save[:,6])
plt.plot(GPS[:,5])
plt.title('东向速度滤波结果')
plt.rcParams['font.sans-serif']=['SimHei']
plt.rcParams['axes.unicode_minus']=False
plt.xlabel('时间')
plt.ylabel('东向速度/(m/s)')
plt.legend(['滤波轨迹','GPS原始数据'])
plt.grid(True)

plt.subplot(312)
plt.plot(data_save[:,7])
plt.plot(GPS[:,6])
plt.title('北向速度滤波结果')
plt.rcParams['font.sans-serif']=['SimHei']
plt.rcParams['axes.unicode_minus']=False
plt.xlabel('时间')
plt.ylabel('北向速度/(m/s)')
plt.legend(['滤波轨迹','GPS原始数据'])
plt.grid(True)

plt.subplot(313)
plt.plot(data_save[:,8])
plt.plot(GPS[:,7])
plt.title('天向速度滤波结果')
plt.rcParams['font.sans-serif']=['SimHei']
plt.rcParams['axes.unicode_minus']=False
plt.xlabel('时间')
plt.ylabel('天向速度/(m/s)')
plt.legend(['滤波轨迹','GPS原始数据'])
plt.grid(True)
plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.2, hspace=0.4)
plt.show()

fig = plt.figure('状态估计阵P阵对角线元素变化曲线图',figsize=(16,9))
plt.subplot(5,3,1)
plt.plot(P_save[:,0])
plt.title('东向失准角')
plt.grid(True)

plt.subplot(5,3,2)
plt.plot(P_save[:,1])
plt.title('北向失准角')
plt.grid(True)

plt.subplot(5,3,3)
plt.plot(P_save[:,2])
plt.title('天向失准角')
plt.grid(True)

plt.subplot(5,3,4)
plt.plot(P_save[:,3])
plt.title('东向速度误差')
plt.grid(True)


plt.subplot(5,3,5)
plt.plot(P_save[:,4])
plt.title('北向速度误差')
plt.grid(True)

plt.subplot(5,3,6)
plt.plot(P_save[:,5])
plt.title('天向速度误差')
plt.grid(True)

plt.subplot(5,3,7)
plt.plot(P_save[:,6])
plt.title('纬度误差')
plt.grid(True)

plt.subplot(5,3,8)
plt.plot(P_save[:,7])
plt.title('经度误差')
plt.grid(True)

plt.subplot(5,3,9)
plt.plot(P_save[:,8])
plt.title('高度误差')
plt.grid(True)

plt.subplot(5,3,10)
plt.plot(P_save[:,9])
plt.title('X轴陀螺漂移')
plt.grid(True)

plt.subplot(5,3,11)
plt.plot(P_save[:,10])
plt.title('Y轴陀螺漂移')
plt.grid(True)

plt.subplot(5,3,12)
plt.plot(P_save[:,11])
plt.title('Z轴陀螺漂移')
plt.grid(True)

plt.subplot(5,3,13)
plt.plot(P_save[:,12])
plt.title('X轴加计零偏')
plt.grid(True)

plt.subplot(5,3,14)
plt.plot(P_save[:,13])
plt.title('Y轴加计零偏')
plt.grid(True)

plt.subplot(5,3,15)
plt.plot(P_save[:,14])
plt.title('Z轴加计零偏')
plt.grid(True)

plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.13, hspace=0.45)
plt.show()

