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

# figure('name', '3D滤波结果');
# plot3(data_save(:,1)*deg,data_save(:,2)*deg,data_save(:,3));
# hold on;
# plot3(GPS(1:i,3)*deg,GPS(1:i,4)*deg,GPS(1:i,5));
# hold on;
# plot3(sins_save(1:i,1)*deg,sins_save(1:i,2)*deg,sins_save(1:i,3));
# title('3D滤波结果');
# xlabel('纬度/°');
# ylabel('经度/°');
# zlabel('高度/h');
# legend('滤波轨迹','GPS原始数据','纯惯性轨迹');
# grid on;