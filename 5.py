import math
import numpy as np

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
# print(GPS.shape)
# print(IMU)

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
g_b = (-np.mean(IMU[0:5000, 5:8], 0) * g).T
# w_ib_b = mean(IMU(1:5000,3:5))'
# A = [g_n';w_ie_n';cross(g_n',w_ie_n')]
# B = [g_b';w_ib_b';cross(g_b',w_ib_b')]
# bCn =A\B
# yaw = atan2(bCn(1,2), bCn(2,2))
# roll = atan2(-bCn(3,1), bCn(3,3))
# pitch = asin(bCn(3,2))

# q = Euler2Quaternion(pitch, roll, yaw)
print(g_b)