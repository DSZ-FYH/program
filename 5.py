import math
import numpy as np
from numpy.linalg import lstsq  # 解超定方程
from numpy.linalg import solve  # 解线性方程
from tqdm import tqdm


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


def Quaternion2bCn(q):
    C = np.zeros((3, 3))
    C[0, 0] = q[0][0] * q[0][0] + q[1][0] * q[1][0] - q[2][0] * q[2][0] - q[3][
        0] * q[3][0]
    C[0, 1] = 2 * (q[1][0] * q[2][0] - q[0][0] * q[3][0])
    C[0, 2] = 2 * (q[1][0] * q[3][0] + q[0][0] * q[2][0])
    C[1, 0] = 2 * (q[1][0] * q[2][0] + q[0][0] * q[3][0])
    C[1, 1] = q[0][0] * q[0][0] - q[1][0] * q[1][0] + q[2][0] * q[2][0] - q[3][
        0] * q[3][0]
    C[1, 2] = 2 * (q[2][0] * q[3][0] - q[0][0] * q[1][0])
    C[2, 0] = 2 * (q[1][0] * q[3][0] - q[0][0] * q[2][0])
    C[2, 1] = 2 * (q[2][0] * q[3][0] + q[0][0] * q[1][0])
    C[2, 2] = q[0][0] * q[0][0] - q[1][0] * q[1][0] - q[2][0] * q[2][0] + q[3][
        0] * q[3][0]
    return C

def Euler2bCn(pitch, roll, yaw):
    C = np.zeros((3, 3))
    C[0, 0] = math.cos(roll) * math.cos(yaw) + math.sin(roll) * math.sin(yaw) * math.sin(pitch)
    C[1, 0] = -math.cos(roll) * math.sin(yaw) + math.sin(roll) * math.cos(yaw) * math.sin(pitch)
    C[2, 0] = -math.sin(roll) * math.cos(pitch)
    C[0, 1] = math.sin(yaw) * math.cos(pitch)
    C[1, 1] = math.cos(yaw) * math.cos(pitch)
    C[2, 1] = math.sin(pitch)
    C[0, 2] = math.sin(roll) * math.cos(yaw) - math.cos(roll) * math.sin(yaw) * math.sin(pitch)
    C[1, 2] = -math.sin(roll) * math.sin(yaw) - math.cos(roll) * math.cos(yaw) * math.sin(pitch)
    C[2, 2] = math.cos(roll) * math.cos(pitch)
    return C

print ("==============参数初始化==============")
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
g_n_transpose = g_n.reshape(1, g_n.shape[0])
w_ie_n = np.array([[0], [Wie * math.cos(lat)], [Wie * math.sin(lat)]])
w_ie_n_transpose = w_ie_n.reshape(1, w_ie_n.shape[0])

g_b_transpose = -np.mean(IMU[0:5000, 5:8], 0) * g
g_b = g_b_transpose.reshape(g_b_transpose.shape[0], 1)
w_ib_b_transpose = np.mean(IMU[0:5000, 2:5], 0)
w_ib_b = w_ib_b_transpose.reshape(w_ib_b_transpose.shape[0], 1)

A = np.array([
    g_n_transpose, w_ie_n_transpose,
    np.cross(g_n_transpose, w_ie_n_transpose)
]).reshape((3, 3))

B = np.array([
    g_b_transpose, w_ib_b_transpose,
    np.cross(g_b_transpose, w_ib_b_transpose)
]).reshape((3, 3))

bCn = solve(A, B)

yaw = math.atan2(bCn[0, 1], bCn[1, 1])
roll = math.atan2(-bCn[2, 0], bCn[2, 2])
pitch = math.asin(bCn[2, 1])

q = Euler2Quaternion(pitch, roll, yaw)

q_sins = q
bCn_sins = bCn
lon_sins = lon
lat_sins = lat
h_sins = h
ve_sins = ve
vn_sins = vn
vu_sins = vu

X_k = np.zeros((15, 1))

array_temp_1 = np.column_stack((np.deg2rad(1) * np.ones(
    (1, 3)), 0.01 * np.ones((1, 3)), (0.1 / Re), (0.1 / Re), 0.15,
                                np.deg2rad(0.1 / 3600) * np.ones(
                                    (1, 3)), (50e-6 * g0) * np.ones(
                                        (1, 3)))).flatten()
P_k = np.dot(np.diag(array_temp_1), np.diag(array_temp_1))

data_save = np.zeros((kalman_end, 9))
X_save = np.zeros((kalman_end, 15))
P_save = np.zeros((kalman_end, 15))
sins_save = np.zeros((kalman_end, 3))

print ("==============初始化完成==============")

print ("开始解算...")

for i in tqdm(range(0, kalman_end)):
    # 计算Rm和Rn
    Rm = Re * (1 - 2 * f + 3 * f * math.sin(lat)**2)
    Rn = Re * (1 + f * math.sin(lat)**2)

    # 去除有害角速度
    w_ie_n = np.array([[0], [Wie * math.cos(lat)], [Wie * math.sin(lat)]])
    w_en_n = np.array([[-vn / (Rm + h)], [ve / (Rn + h)],
                       [ve / (Rn + h) * math.tan(lat)]])
    w_in_n = w_ie_n + w_en_n
    w_in_b = np.matmul(bCn.T, w_in_n)
    w_ib_b = IMU[i, 2:5].reshape(IMU[i, 2:5].shape[0], 1) - gyro_bias
    w_nb_b = w_ib_b - w_in_b

    # 四元数更新
    dTheta = w_nb_b * kalman_gap * IMU_dT
    dTheta_norm = np.linalg.norm(dTheta)

    q1 = np.array([
        np.float64(math.cos(dTheta_norm / 2)),
        float(-dTheta[0][0] / dTheta_norm * math.sin(dTheta_norm / 2)),
        float(-dTheta[1][0] / dTheta_norm * math.sin(dTheta_norm / 2)),
        float(-dTheta[2][0] / dTheta_norm * math.sin(dTheta_norm / 2))
    ])
    q2 = np.array([
        float(dTheta[0][0] / dTheta_norm * math.sin(dTheta_norm / 2)),
        np.float64(math.cos(dTheta_norm / 2)),
        float(dTheta[2][0] / dTheta_norm * math.sin(dTheta_norm / 2)),
        float(-dTheta[1][0] / dTheta_norm * math.sin(dTheta_norm / 2))
    ])
    q3 = np.array([
        float(dTheta[1][0] / dTheta_norm * math.sin(dTheta_norm / 2)),
        float(-dTheta[2][0] / dTheta_norm * math.sin(dTheta_norm / 2)),
        np.float64(math.cos(dTheta_norm / 2)),
        float(dTheta[0][0] / dTheta_norm * math.sin(dTheta_norm / 2))
    ])
    q4 = np.array([
        float(dTheta[2][0] / dTheta_norm * math.sin(dTheta_norm / 2)),
        float(dTheta[1][0] / dTheta_norm * math.sin(dTheta_norm / 2)),
        float(-dTheta[0][0] / dTheta_norm * math.sin(dTheta_norm / 2)),
        np.float64(math.cos(dTheta_norm / 2))
    ])
    q = np.matmul(np.array([q1, q2, q3, q4]), q)
    q = q / np.linalg.norm(q)  # 单位化四元数
    bCn = Quaternion2bCn(q)  # 更新Cbn阵

    # 重力加速度更新
    g = g0 * (1 + 0.00193185138639 * math.sin(lat)**2) / math.sqrt(
        1 - 0.00669437999013 * math.sin(lat)**2) * Re**2 / (Re + h)**2

    # 速度更新
    fb = (IMU[i, 5:8] * 9.8).reshape((IMU[i, 5:8].shape[0], 1)) - acc_bias
    fn = np.matmul(bCn, fb)

    w_ie_n = w_ie_n.reshape((-1, 1))
    w_en_n = w_en_n.reshape((-1, 1))

    # 比力方程
    dv = np.array(fn - np.cross((2 * w_ie_n + w_en_n).reshape(
        (1, -1)).flatten(), np.array([ve, vn, vu]).reshape((1,-1)).flatten()).reshape(
            (-1, 1)) + np.array([0, 0, -g]).reshape(
                (np.array([0, 0, -g]).shape[0], 1))) * kalman_gap * IMU_dT

    # print(dv)

    ve = ve + dv[0]
    vn = vn + dv[1]
    vu = vu + dv[2]

    # 位置更新
    lat = lat + vn * kalman_gap * IMU_dT / (Rm + h)
    lon = lon + ve * kalman_gap * IMU_dT / ((Rn + h) * math.cos(lat))
    h = h + vu * kalman_gap * IMU_dT

    # 纯惯性捷联解算
    w_ie_n_sins = np.array(
        [0, Wie * math.cos(lat_sins), Wie * math.sin(lat_sins)]).reshape(
            (-1, 1))
    w_en_n_sins = np.array([
        -vn_sins / (Rm + h_sins), ve_sins / (Rn + h_sins),
        ve_sins / (Rn + h_sins) * math.tan(lat_sins)
    ]).reshape((-1, 1))
    w_in_n_sins = w_ie_n_sins + w_en_n_sins
    w_in_b_sins = np.matmul(bCn_sins.transpose(), w_in_n_sins)
    w_ib_b_sins = IMU[i, 2:5].reshape((-1, 1))
    w_nb_b_sins = w_ib_b_sins - w_in_b_sins

    dTheta = w_nb_b_sins*kalman_gap*IMU_dT
    dTheta_norm = np.linalg.norm(dTheta)

    q_sins1 = np.array([
        np.float64(math.cos(dTheta_norm / 2)),
        float(-dTheta[0][0] / dTheta_norm * math.sin(dTheta_norm / 2)),
        float(-dTheta[1][0] / dTheta_norm * math.sin(dTheta_norm / 2)),
        float(-dTheta[2][0] / dTheta_norm * math.sin(dTheta_norm / 2))
    ])
    q_sins2 = np.array([
        float(dTheta[0][0] / dTheta_norm * math.sin(dTheta_norm / 2)),
        np.float64(math.cos(dTheta_norm / 2)),
        float(dTheta[2][0] / dTheta_norm * math.sin(dTheta_norm / 2)),
        float(-dTheta[1][0] / dTheta_norm * math.sin(dTheta_norm / 2))
    ])
    q_sins3 = np.array([
        float(dTheta[1][0] / dTheta_norm * math.sin(dTheta_norm / 2)),
        float(-dTheta[2][0] / dTheta_norm * math.sin(dTheta_norm / 2)),
        np.float64(math.cos(dTheta_norm / 2)),
        float(dTheta[0][0] / dTheta_norm * math.sin(dTheta_norm / 2))
    ])
    q_sins4 = np.array([
        float(dTheta[2][0] / dTheta_norm * math.sin(dTheta_norm / 2)),
        float(dTheta[1][0] / dTheta_norm * math.sin(dTheta_norm / 2)),
        float(-dTheta[0][0] / dTheta_norm * math.sin(dTheta_norm / 2)),
        np.float64(math.cos(dTheta_norm / 2))
    ])
    q_sins = np.matmul(np.array([q_sins1, q_sins2, q_sins3, q_sins4]), q_sins)

    q_sins = q_sins/np.linalg.norm(q_sins)      # 单位化四元数
    bCn_sins = Quaternion2bCn(q_sins)           # 更新Cbn

    fb_sins = IMU[i,5:8].reshape((-1,1))
    fn_sins = np.matmul(bCn_sins,fb)

    w_ie_n_sins = w_ie_n_sins.reshape((-1, 1))
    w_en_n_sins = w_en_n_sins.reshape((-1, 1))

    dv = np.array(fn_sins - np.cross((2 * w_ie_n_sins + w_en_n_sins).reshape(
        (1, -1)).flatten(), np.array([ve_sins, vn_sins, vu_sins]).reshape((1,-1)).flatten()).reshape(
            (-1, 1)) + np.array([0, 0, -g]).reshape(
                (np.array([0, 0, -g]).shape[0], 1))) * kalman_gap * IMU_dT # 比力方程 
    ve_sins = ve_sins + dv[0]
    vn_sins = vn_sins + dv[1]
    vu_sins = vu_sins + dv[2]
    
    lat_sins = lat_sins + vn*kalman_gap*IMU_dT/(Rm+h)
    lon_sins = lon_sins + ve*kalman_gap*IMU_dT/((Rn+h)*math.cos(lat))
    h_sins = h_sins + vu*kalman_gap*IMU_dT

    # 纯惯性数据保存
    lat_sins = lat_sins[0]
    lon_sins = lon_sins[0]
    h_sins = h_sins[0]
    sins_save[i,0:3] = [lat_sins,lon_sins,h_sins]

    # 卡尔曼模型更新
    F = np.zeros((15,15))
    F[0,1] = Wie*math.sin(lat)+ ve/(Rn+h)*math.tan(lat)
    F[0,2] = -(Wie*math.cos(lat)+ve/(Rn+h))
    F[0,4] = -1/(Rm+h)
    F[1,0] = -(Wie*math.sin(lat)+ve/(Rn+h)*math.tan(lat))
    F[1,2] = -vn/(Rm+h)
    F[1,3] = 1/(Rn+h)
    F[1,6] = -Wie*math.sin(lat)
    F[2,0] = Wie*math.cos(lat)+ve/(Rn+h)
    F[2,1] = vn/(Rm+h)
    F[2,3] = 1/(Rn+h)*math.tan(lat)
    F[2,6] = Wie*math.cos(lat)+ve/(Rn+h)*1/math.cos(lat)**2
    F[3,1] = -fn[2]     # f_u天向比力
    F[3,2] = fn[1]      # f_n北向比力
    F[3,3] = (vn*math.tan(lat)-vu)/(Rn+h)
    F[3,4] = 2*Wie*math.sin(lat)+ve/(Rn+h)*math.tan(lat)
    F[3,5] = -(2*Wie*math.cos(lat)+ve/(Rn+h))
    F[3,6] = 2*Wie*math.cos(lat)*vn+ve*vn/(Rn+h)*1/math.cos(lat)**2+2*Wie*math.sin(lat)*vu
    F[4,0] = fn[2]      # f_u天向比力
    F[4,2] = -fn[0]     # f_e东向比力
    F[4,3] = -(2*Wie*math.sin(lat)+ve/(Rn+h)*math.tan(lat))
    F[4,4] = -vu/(Rm+h)
    F[4,5] = -vn/(Rm+h)
    F[4,6] = -(2*Wie*math.cos(lat)+ve/(Rn+h)*1/math.cos(lat)**2)*ve
    F[5,0] = -fn[1]     # f_n北向比力
    F[5,1] = fn[0]      # f_e东向比力
    F[5,3] = 2*(Wie*math.cos(lat)+ve/(Rn+h))
    F[5,4] = 2*vn/(Rm+h)
    F[5,6] = -2*ve*Wie*math.sin(lat)
    F[6,4] = 1/(Rm+h)
    F[7,3] = 1/math.cos(lat)/(Rn+h)
    F[7,6]= ve/(Rn+h)*1/math.cos(lat)*math.tan(lat)
    F[8,5]  = 1
    F[0:3,9:12] = bCn
    F[3:6,12:15] = bCn

    # 状态转移矩阵F离散化
    F = np.eye(15)+F*IMU_dT+np.linalg.matrix_power(F,2)*IMU_dT**2/2+np.linalg.matrix_power(F,3)*IMU_dT**3/6

    G = np.zeros((15,6))
    G[0:3,0:3] = bCn
    G[3:6,3:6] = bCn

    H = np.zeros((6,15))
    H[0:3,3:6] = np.eye(3)
    H[3:6,6:9] = np.eye(3)
    # H[0:6,3:9] = np.eye(6)

    Q = np.linalg.matrix_power(np.diag(np.array([np.ones(3)*np.deg2rad(0.1/3600),np.ones(3)*50e-6*g]).reshape(1,-1).flatten()),2)

    R_temp = np.array([np.ones(1)*0.1/(Rn*math.cos(lat)),np.ones(1)*0.1/Rm,np.ones(1)*0.15]).reshape(1,-1).flatten()
    R = np.linalg.matrix_power(np.diag(np.array([np.ones(3)*0.01,R_temp]).reshape(1,-1).flatten()),2)

    Z =np.append(np.array([ve,vn,vu]).reshape(1,-1).flatten()-GPS[i,5:8],np.array([lat,lon,h]).reshape(1,-1).flatten()-GPS[i,2:5]).reshape(-1,1)

    X_k_k_1 = np.matmul(F,X_k)
    P_k_k_1 = np.matmul(np.matmul(F,P_k),F.T)+np.matmul(np.matmul(G,Q),G.T)
    K_k = np.matmul(np.matmul(P_k_k_1,H.T),np.linalg.inv(np.matmul(np.matmul(H,P_k_k_1),H.T)+R))
    X_k = X_k_k_1 + np.matmul(K_k,Z-np.matmul(H,X_k_k_1))
    P_k = np.matmul(np.eye(15)-np.matmul(K_k,H),P_k_k_1)

    # 状态量（误差）保存
    X_save[i,:] = X_k.T

    # P阵对角线保存
    P_save[i,:] = np.diag(P_k).reshape(1,-1)

    # 姿态修正
    bCn = np.matmul(Euler2bCn(X_k[0],X_k[1],X_k[2]),bCn)
    yaw = np.rad2deg(math.atan2(bCn[0,1], bCn[1,1]))
    roll = np.rad2deg(math.atan2(-bCn[2,0], bCn[2,2]))
    pitch = np.rad2deg(math.asin(bCn[2,1]))

    # 速度修正
    ve = ve[0] - X_k[3][0]
    vn = vn[0] - X_k[4][0]
    vu = vu[0] - X_k[5][0]

    # 位置修正
    lat = lat[0] - X_k[6][0]
    lon = lon[0] - X_k[7][0]
    h = h[0] - X_k[8][0]

    gyro_bias = X_k[9:12]
    acc_bias = X_k[12:15]

    # 误差清零
    X_k[0:9] = np.zeros((9,1))
    
    # 滤波数据保存
    data_save[i,:] = np.array([lat,lon,h,pitch,roll,yaw,ve,vn,vu])
    

np.savetxt("data_save.txt",data_save,delimiter=" ")
np.savetxt("GPS.txt",GPS,delimiter=" ")
np.savetxt("sins_save.txt",sins_save,delimiter=" ")
np.savetxt("P_save.txt",P_save,delimiter=" ")

print ("解算完成！")
print ("开始绘图")