import math
import numpy as np
from numpy.linalg import lstsq  # 解超定方程
from numpy.linalg import solve  # 解线性方程


def Euler2Quaternion(pitch, roll, yaw):
    q = np.array([
        math.cos(yaw / 2) * math.cos(pitch / 2) * math.cos(roll / 2) -
        math.sin(yaw / 2) * math.sin(pitch / 2) * math.sin(roll / 2),
        math.cos(yaw / 2) * math.sin(pitch / 2) * math.cos(roll / 2) -
        math.sin(yaw / 2) * math.cos(pitch / 2) * math.sin(roll / 2),
        math.cos(yaw / 2) * math.cos(pitch / 2) * math.sin(roll / 2) +
        math.sin(yaw / 2) * math.sin(pitch / 2) * math.cos(roll / 2),
        math.cos(yaw / 2) * math.sin(pitch / 2) * math.sin(roll / 2) +
        math.sin(yaw / 2) * math.cos(pitch / 2) * math.cos(roll / 2)
    ]).reshape((4, 1))
    return q


# 参数设置
IMU_frequency = 200
GPS_frequency = 1
IMU_dT = 1 / IMU_frequency
GPS_dT = 1 / GPS_frequency
rad = math.pi / 180
deg = 180 / math.pi
g0 = 9.7803267714  # 重力加速度
Re = 6378137  # 地球半长轴
Rp = 6356752
f = (Re - Rp) / Re
e = 1 / 298.257
Rm = Re
Rn = Re
Wie = 15 / 3600 * rad  # 地球自转角速度（单位：弧度/h）

# 导入数据
GPS = np.loadtxt("实验5-GPS-INS组合跑车实验\GPS插值数据.dat")
IMU = np.loadtxt("实验5-GPS-INS组合跑车实验\惯导数据.dat")
GPS_data_length = np.size(GPS, 0)
IMU_data_length = np.size(IMU, 0)

# 滤波起始停止长度设置
# matlab 长度设置可能不同，循环可能有问题
kalman_start = 0
kalman_end = IMU_data_length
kalman_gap = 1

# 数据单位转换
GPS[:, 2:4] = GPS[:, 2:4] * rad  # 经纬度
IMU[:, 2:5] = IMU[:, 2:5] / 3600 * rad  # 角速度

# 初值设置
lon = GPS[kalman_start][3]
lat = GPS[kalman_start][2]
h = GPS[kalman_start][4]
ve = GPS[kalman_start][5]
vn = GPS[kalman_start][6]
vu = GPS[kalman_start][7]
gyro_bias = np.zeros((3, 1), dtype=np.float32)
acc_bias = np.zeros((3, 1), dtype=np.float32)

# 粗对准
g = g0 * (1 + 0.00193185138639 * math.sin(lat)**2) / math.sqrt(
    1 - 0.00669437999013 * math.sin(lat)**2) * Re**2 / (Re + h)**2
g_n = np.array([[0], [0], [-g]])
w_ie_n = np.array([[0], [Wie * math.cos(lat)], [Wie * math.sin(lat)]])
g_b = -np.mean(IMU[0:5000, 5:8], 0) * g
g_b = g_b.reshape(g_b.shape[0], 1)
w_ib_b = np.mean(IMU[0:5000, 2:5], 0)
w_ib_b = w_ib_b.reshape(w_ib_b.shape[0], 1)

g_n_transpose = g_n.reshape(1, g_n.shape[0])
w_ie_n_transpose = w_ie_n.reshape(1, w_ie_n.shape[0])

A = np.array([
    g_n_transpose, w_ie_n_transpose,
    np.cross(g_n_transpose, w_ie_n_transpose)
]).reshape((3, 3))

g_b_transpose = g_b.reshape(1, g_b.shape[0])
w_ib_b_transpose = w_ib_b.reshape(1, w_ib_b.shape[0])

B = np.array([
    g_b_transpose, w_ib_b_transpose,
    np.cross(g_b_transpose, w_ib_b_transpose)
]).reshape((3, 3))

bCn = solve(A, B)

yaw = math.atan2(bCn[0, 1], bCn[1, 1])
roll = math.atan2(-bCn[2, 0], bCn[2, 2])
pitch = math.asin(bCn[2, 1])

q = Euler2Quaternion(pitch, roll, yaw)
