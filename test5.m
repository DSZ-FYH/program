close all;
clear;
clc;

format long g;
format compact;

% 参数设置
IMU_frequency = 200;
GPS_frequency = 1;
IMU_dT = 1/IMU_frequency;
GPS_dT = 1/GPS_frequency;
rad = pi/180;
deg = 180/pi;
g0 = 9.7803267714;
Re = 6378137; %地球半长轴
Rp = 6356752;
f = (Re-Rp)/Re;
e = 1/298.257;
Rm = Re;
Rn = Re;
Wie = 15/3600*rad; %地球自转角速度（单位：弧度/h）

% 导入数据
GPS = importdata('.\实验5-GPS-INS组合跑车实验\GPS数据.dat');
IMU = importdata('.\实验5-GPS-INS组合跑车实验\惯导数据.dat');
GPS_data_length = length(GPS);
IMU_data_length = length(IMU);

% GPS数据差值
GPS = interp1(1:GPS_data_length, GPS, 1:0.005:GPS_data_length, 'spline');
GPS(1,:) = [];

% 滤波起始停止长度设置
kalman_start = 1;
kalman_end = IMU_data_length;
kalman_gap = 1;

% 数据单位转换
GPS(:,3:4) = GPS(:,3:4)*rad; %经纬度
IMU(:,3:5) = IMU(:,3:5)/3600*rad; %角速度

% 初值设置
lon = GPS(kalman_start,4);
lat = GPS(kalman_start,3);
h = GPS(kalman_start,5);
ve = GPS(kalman_start,6);
vn = GPS(kalman_start,7);
vu = GPS(kalman_start,8);
gyro_bias = [0;0;0];
acc_bias = [0;0;0];

i=0

% 粗对准
g = g0*(1+0.00193185138639*sin(lat)^2)/sqrt(1-0.00669437999013*sin(lat)^2)*Re^2/(Re+h)^2;
g_n = [0;0;-g];
w_ie_n = [0;Wie*cos(lat);Wie*sin(lat)];
g_b = -mean(IMU(1:5000,6:8))'*g;
w_ib_b = mean(IMU(1:5000,3:5))';
A = [g_n';w_ie_n';cross(g_n',w_ie_n')];
B = [g_b';w_ib_b';cross(g_b',w_ib_b')];
bCn =A\B;
yaw = atan2(bCn(1,2), bCn(2,2));
roll = atan2(-bCn(3,1), bCn(3,3));
pitch = asin(bCn(3,2));

q = Euler2Quaternion(pitch, roll, yaw);

i=0;

q_sins = q;
bCn_sins = bCn;
lon_sins = lon;
lat_sins = lat;
h_sins = h;
ve_sins = ve;
vn_sins = vn;
vu_sins = vu;

X_k = zeros(15,1);
P_k = diag([deg2rad(1)*ones(1,3),0.01*ones(1,3),(0.1/Re),(0.1/Re),0.15,deg2rad(0.1/3600)*ones(1,3),(50e-6*g0)*ones(1,3)])^2;

data_save = zeros(kalman_end,9);
X_save = zeros(kalman_end,15);
P_save = zeros(kalman_end,15);
sins_save = zeros(kalman_end,3);

for i=kalman_start:kalman_gap:kalman_end
    
    % i
    
    % 计算Rm和Rn
    Rm = Re*(1-2*f+3*f*sin(lat)^2);
    Rn = Re*(1+f*sin(lat)^2);
    
    % 去除有害角速度
    w_ie_n = [0;Wie*cos(lat);Wie*sin(lat)];
    w_en_n = [-vn/(Rm+h);ve/(Rn+h);ve/(Rn+h)*tan(lat)];
    w_in_n = w_ie_n + w_en_n;
    w_in_b = bCn'*w_in_n;
    w_ib_b = IMU(i,3:5)' - gyro_bias;
    w_nb_b = w_ib_b - w_in_b;
    
    % 四元数更新
    dTheta = w_nb_b*kalman_gap*IMU_dT;
    dTheta_norm = norm(dTheta);
    q = [cos(dTheta_norm/2) -dTheta(1)/dTheta_norm*sin(dTheta_norm/2) -dTheta(2)/dTheta_norm*sin(dTheta_norm/2) -dTheta(3)/dTheta_norm*sin(dTheta_norm/2);
         dTheta(1)/dTheta_norm*sin(dTheta_norm/2) cos(dTheta_norm/2) dTheta(3)/dTheta_norm*sin(dTheta_norm/2) -dTheta(2)/dTheta_norm*sin(dTheta_norm/2);
         dTheta(2)/dTheta_norm*sin(dTheta_norm/2) -dTheta(3)/dTheta_norm*sin(dTheta_norm/2) cos(dTheta_norm/2) dTheta(1)/dTheta_norm*sin(dTheta_norm/2);
         dTheta(3)/dTheta_norm*sin(dTheta_norm/2) dTheta(2)/dTheta_norm*sin(dTheta_norm/2) -dTheta(1)/dTheta_norm*sin(dTheta_norm/2) cos(dTheta_norm/2)]*q;
%     q = 0.5*dT*[0 -w_nb_b(1) -w_nb_b(2) -w_nb_b(3);w_nb_b(1) 0 w_nb_b(3) -w_nb_b(2);w_nb_b(2) -w_nb_b(3) 0 w_nb_b(1);w_nb_b(3) w_nb_b(2) -w_nb_b(1) 0]*q+q;
	q = q/norm(q); %单位化四元数
    bCn = Quaternion2bCn(q); %更新Cbn阵
    
    % 重力加速度更新
    g = g0*(1+0.00193185138639*sin(lat)^2)/sqrt(1-0.00669437999013*sin(lat)^2)*Re^2/(Re+h)^2;
    
    % 速度更新
    fb = IMU(i,6:8)'*9.8 - acc_bias;
    fn = bCn*fb;
    dv = (fn - cross((2*w_ie_n+w_en_n),[ve;vn;vu]) + [0;0;-g])*kalman_gap*IMU_dT; %比力方程
    ve = ve + dv(1);
    vn = vn + dv(2);
    vu = vu + dv(3);
    
    % 位置更新
    lat = lat + vn*kalman_gap*IMU_dT/(Rm+h);
    lon = lon + ve*kalman_gap*IMU_dT/((Rn+h)*cos(lat));
    h = h + vu*kalman_gap*IMU_dT;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 纯惯性捷联解算
    w_ie_n_sins = [0;Wie*cos(lat_sins);Wie*sin(lat_sins)];
    w_en_n_sins = [-vn_sins/(Rm+h_sins);ve_sins/(Rn+h_sins);ve_sins/(Rn+h_sins)*tan(lat_sins)];
    w_in_n_sins = w_ie_n_sins + w_en_n_sins;
    w_in_b_sins = bCn_sins'*w_in_n_sins;
    w_ib_b_sins = IMU(i,3:5)';
    w_nb_b_sins = w_ib_b_sins - w_in_b_sins;
    
    dTheta = w_nb_b_sins*kalman_gap*IMU_dT;
    dTheta_norm = norm(dTheta);
    q_sins = [cos(dTheta_norm/2) -dTheta(1)/dTheta_norm*sin(dTheta_norm/2) -dTheta(2)/dTheta_norm*sin(dTheta_norm/2) -dTheta(3)/dTheta_norm*sin(dTheta_norm/2);
         dTheta(1)/dTheta_norm*sin(dTheta_norm/2) cos(dTheta_norm/2) dTheta(3)/dTheta_norm*sin(dTheta_norm/2) -dTheta(2)/dTheta_norm*sin(dTheta_norm/2);
         dTheta(2)/dTheta_norm*sin(dTheta_norm/2) -dTheta(3)/dTheta_norm*sin(dTheta_norm/2) cos(dTheta_norm/2) dTheta(1)/dTheta_norm*sin(dTheta_norm/2);
         dTheta(3)/dTheta_norm*sin(dTheta_norm/2) dTheta(2)/dTheta_norm*sin(dTheta_norm/2) -dTheta(1)/dTheta_norm*sin(dTheta_norm/2) cos(dTheta_norm/2)]*q_sins;
%     q = 0.5*dT*[0 -w_nb_b(1) -w_nb_b(2) -w_nb_b(3);w_nb_b(1) 0 w_nb_b(3) -w_nb_b(2);w_nb_b(2) -w_nb_b(3) 0 w_nb_b(1);w_nb_b(3) w_nb_b(2) -w_nb_b(1) 0]*q+q;
	q_sins = q_sins/norm(q_sins); %单位化四元数
    bCn_sins = Quaternion2bCn(q_sins); %更新Cbn
    
    fb_sins = IMU(i,6:8)';
    fn_sins = bCn_sins*fb;
    dv = (fn_sins - cross((2*w_ie_n_sins+w_en_n_sins),[ve_sins;vn_sins;vu_sins]) + [0;0;-g])*kalman_gap*IMU_dT; %比力方程
    ve_sins = ve_sins + dv(1);
    vn_sins = vn_sins + dv(2);
    vu_sins = vu_sins + dv(3);
    
    lat_sins = lat_sins + vn*kalman_gap*IMU_dT/(Rm+h);
    lon_sins = lon_sins + ve*kalman_gap*IMU_dT/((Rn+h)*cos(lat));
    h_sins = h_sins + vu*kalman_gap*IMU_dT;
    
    % 纯惯性数据保存
    sins_save(i,:) = [lat_sins lon_sins h_sins];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    % 卡尔曼模型更新
    F = zeros(15);
    F(1,2) = Wie*sin(lat)+ ve/(Rn+h)*tan(lat);
    F(1,3) = -(Wie*cos(lat)+ve/(Rn+h));
    F(1,5) = -1/(Rm+h);
    F(2,1) = -Wie*sin(lat)-ve/(Rn+h)*tan(lat);
    F(2,3) = -vn/(Rm+h);
    F(2,4) = 1/(Rn+h);
    F(2,7) = -Wie*sin(lat);
    F(3,1) = Wie*cos(lat)+ve/(Rn+h);
    F(3,2) = vn/(Rm+h);
    F(3,4) = 1/(Rn+h)*tan(lat);
    F(3,7) = Wie*cos(lat)+ve/(Rn+h)*sec(lat)^2;
    F(4,2) = -fn(3); %f_u天向比力
    F(4,3) = fn(2); %f_n北向比力
    F(4,4) = (vn*tan(lat)-vu)/(Rn+h);
    F(4,5) = 2*Wie*sin(lat)+ve/(Rn+h)*tan(lat);
    F(4,6) = -(2*Wie*cos(lat)+ve/(Rn+h));
    F(4,7) = 2*Wie*cos(lat)*vn+ve*vn/(Rn+h)*sec(lat)^2+2*Wie*sin(lat)*vu;
    F(5,1) = fn(3); %f_u天向比力
    F(5,3) = -fn(1); %f_e东向比力
    F(5,4) = -(2*Wie*sin(lat)+ve/(Rn+h)*tan(lat));
    F(5,5) = -vu/(Rm+h);
    F(5,6) = -vn/(Rm+h);
    F(5,7) = -(2*Wie*cos(lat)+ve/(Rn+h)*sec(lat)^2)*ve;
    F(6,1) = -fn(2); %f_n北向比力
    F(6,2) = fn(1); %f_e东向比力
    F(6,4) = 2*(Wie*cos(lat)+ve/(Rn+h));
    F(6,5) = 2*vn/(Rm+h);
    F(6,7) = -2*ve*Wie*sin(lat);
    F(7,5) = 1/(Rm+h);
    F(8,4) = sec(lat)/(Rn+h);
    F(8,7) = ve/(Rn+h)*sec(lat)*tan(lat);
    F(9,6) = 1;
    F(1:3,10:12) = bCn;
    F(4:6,13:15) = bCn;
    % 状态转移矩阵F离散化
    F = eye(15)+F*IMU_dT+F^2*IMU_dT^2/factorial(2)+F^3*IMU_dT^3/factorial(3);

    G = zeros(15,6);
    G(1:3,1:3) = bCn;
    G(4:6,4:6) = bCn;

    H = zeros(6,15);
    H(1:3,4:6) = eye(3);
    H(4:6,7:9) = eye(3);
%         H(4:6,7:9) = diag([Rm+h,(Rn+h)*cos(lat),1]);

    Q = diag([ones(3,1)*deg2rad(0.1/3600);ones(3,1)*(50e-6*g)])^2;

    R = diag([0.01*ones(3,1);0.1/(Rn*cos(lat));0.1/Rm;0.15])^2;

    Z = [[ve;vn;vu]-GPS(i,6:8)'; [lat;lon;h]-GPS(i,3:5)'];

    X_k_k_1 = F*X_k;
    P_k_k_1 = F*P_k*F' + G*Q*G';
    K_k = P_k_k_1*H'*(H*P_k_k_1*H'+R)^-1;
    X_k = X_k_k_1 + K_k*(Z-H*X_k_k_1);
    P_k = (eye(15)-K_k*H)*P_k_k_1;

    % 状态量（误差）保存
    X_save(i,:) = X_k';

    % P阵对角线保存
    P_save(i,:) = diag(P_k)';

    % 姿态修正
%     bCn = [1,-X_k(3),X_k(2);X_k(3),1,-X_k(1);-X_k(2),X_k(1),1]*bCn;
    bCn = Euler2bCn(X_k(1),X_k(2),X_k(3))*bCn;
    yaw = atan2d(bCn(1,2), bCn(2,2));
    roll = atan2d(-bCn(3,1), bCn(3,3));
    pitch = asind(bCn(3,2));

    % 速度修正
    ve = ve - X_k(4);
    vn = vn - X_k(5);
    vu = vu - X_k(6);

    % 位置修正
    lat = lat - X_k(7);
    lon = lon - X_k(8);
    h = h - X_k(9);

    gyro_bias = X_k(10:12);
    acc_bias = X_k(13:15);

    % 误差清零
    X_k(1:9) = zeros(9,1);
    
    % 滤波数据保存
    data_save(i,:) = [lat,lon,h,pitch,roll,yaw,ve,vn,vu];
end

figure('name', '3D滤波结果');
plot3(data_save(:,1)*deg,data_save(:,2)*deg,data_save(:,3));
hold on;
plot3(GPS(1:i,3)*deg,GPS(1:i,4)*deg,GPS(1:i,5));
hold on;
plot3(sins_save(1:i,1)*deg,sins_save(1:i,2)*deg,sins_save(1:i,3));
title('3D滤波结果');
xlabel('纬度/°');
ylabel('经度/°');
zlabel('高度/h');
legend('滤波轨迹','GPS原始数据','纯惯性轨迹');
grid on;

figure('name', '2D滤波结果');
plot(data_save(:,1)*deg,data_save(:,2)*deg);
hold on;
plot(GPS(1:i,3)*deg,GPS(1:i,4)*deg);
hold on;
plot(sins_save(1:i,1)*deg,sins_save(1:i,2)*deg);
title('2D滤波结果');
xlabel('纬度/°');
ylabel('经度/°');
legend('滤波轨迹','GPS原始数据','纯惯性轨迹');
grid on;

figure('name', '位置滤波结果');
subplot(3,1,1);
plot(data_save(:,1)*deg);
hold on;
plot(GPS(1:i,3)*deg);
title('纬度滤波结果');
xlabel('时间');
ylabel('纬度/°');
legend('滤波轨迹','GPS原始数据');
grid on;

subplot(3,1,2);
plot(data_save(:,2)*deg);
hold on;
plot(GPS(1:i,4)*deg);
title('经度滤波结果');
xlabel('时间');
ylabel('经度/°');
legend('滤波轨迹','GPS原始数据');
grid on;

subplot(3,1,3);
plot(data_save(:,3));
hold on;
plot(GPS(1:i,5));
title('高度滤波结果');
xlabel('时间');
ylabel('高度/h');
legend('滤波轨迹','GPS原始数据');
grid on;

figure('name', '姿态滤波结果');
subplot(3,1,1);
plot(data_save(:,4));
title('俯仰滤波结果');
xlabel('时间');
ylabel('俯仰/°');
grid on;

subplot(3,1,2);
plot(data_save(:,5));
title('横滚滤波结果');
xlabel('时间');
ylabel('横滚/°');
grid on;

subplot(3,1,3);
plot(data_save(:,6));
title('偏航滤波结果');
xlabel('时间');
ylabel('偏航/°');
grid on;

figure('name', '速度滤波结果');
subplot(3,1,1);
plot(data_save(:,7));
hold on;
plot(GPS(1:i,6));
title('东向速度滤波结果');
xlabel('时间');
ylabel('东向速度/(m/s)');
legend('滤波轨迹','GPS原始数据');
grid on;

subplot(3,1,2);
plot(data_save(:,8));
hold on;
plot(GPS(1:i,7));
title('北向速度滤波结果');
xlabel('时间');
ylabel('北向速度/(m/s)');
legend('滤波轨迹','GPS原始数据');
grid on;

subplot(3,1,3);
plot(data_save(:,9));
hold on;
plot(GPS(1:i,8));
title('天向速度滤波结果');
xlabel('时间');
ylabel('天向速度/(m/s)');
legend('滤波轨迹','GPS原始数据');
grid on;

figure('name','状态估计阵P阵对角线元素变化曲线图');
% suptitle('状态估计阵P阵对角线元素变化曲线图');
subplot(5,3,1);
plot(P_save(:,1));
title('东向失准角');
grid on;

subplot(5,3,2);
plot(P_save(:,2));
title('北向失准角');
grid on;

subplot(5,3,3);
plot(P_save(:,3));
title('天向失准角');
grid on;

subplot(5,3,4);
plot(P_save(:,4));
title('东向速度误差');
grid on;

subplot(5,3,5);
plot(P_save(:,5));
title('北向速度误差');
grid on;

subplot(5,3,6);
plot(P_save(:,6));
title('天向速度误差');
grid on;

subplot(5,3,7);
plot(P_save(:,7));
title('纬度误差');
grid on;

subplot(5,3,8);
plot(P_save(:,8));
title('经度误差');
grid on;

subplot(5,3,9);
plot(P_save(:,9));
title('高度误差');
grid on;

subplot(5,3,10);
plot(P_save(:,10));
title('X轴陀螺漂移');
grid on;

subplot(5,3,11);
plot(P_save(:,11));
title('Y轴陀螺漂移');
grid on;

subplot(5,3,12);
plot(P_save(:,12));
title('Z轴陀螺漂移');
grid on;

subplot(5,3,13);
plot(P_save(:,13));
title('X轴加计零偏');
grid on;

subplot(5,3,14);
plot(P_save(:,14));
title('Y轴加计零偏');
grid on;

subplot(5,3,15);
plot(P_save(:,15));
title('Z轴加计零偏');
grid on;

sgtitle('状态估计阵P阵对角线元素变化曲线图');

function C = Euler2nCb(pitch, roll, yaw)
    C = zeros(3);
    C(1,1) = cos(roll)*cos(yaw) + sin(roll)*sin(yaw)*sin(pitch);
    C(1,2) = -cos(roll)*sin(yaw) + sin(roll)*cos(yaw)*sin(pitch);
    C(1,3) = -sin(roll)*cos(pitch);
    C(2,1) = sin(yaw)*cos(pitch);
    C(2,2) = cos(yaw)*cos(pitch);
    C(2,3) = sin(pitch);
    C(3,1) = sin(roll)*cos(yaw) - cos(roll)*sin(yaw)*sin(pitch);
    C(3,2) = -sin(roll)*sin(yaw) - cos(roll)*cos(yaw)*sin(pitch);
    C(3,3) = cos(roll)*cos(pitch);
end

function C = Euler2bCn(pitch, roll, yaw)
    C = zeros(3);
    C(1,1) = cos(roll)*cos(yaw) + sin(roll)*sin(yaw)*sin(pitch);
    C(2,1) = -cos(roll)*sin(yaw) + sin(roll)*cos(yaw)*sin(pitch);
    C(3,1) = -sin(roll)*cos(pitch);
    C(1,2) = sin(yaw)*cos(pitch);
    C(2,2) = cos(yaw)*cos(pitch);
    C(3,2) = sin(pitch);
    C(1,3) = sin(roll)*cos(yaw) - cos(roll)*sin(yaw)*sin(pitch);
    C(2,3) = -sin(roll)*sin(yaw) - cos(roll)*cos(yaw)*sin(pitch);
    C(3,3) = cos(roll)*cos(pitch);
end

function C = Quaternion2bCn(q)
    C = zeros(3);
    C(1,1) = q(1)^2+q(2)^2-q(3)^2-q(4)^2;
    C(1,2) = 2*(q(2)*q(3)-q(1)*q(4));
    C(1,3) = 2*(q(2)*q(4)+q(1)*q(3));
    C(2,1) = 2*(q(2)*q(3)+q(1)*q(4));
    C(2,2) = q(1)^2-q(2)^2+q(3)^2-q(4)^2;
    C(2,3) = 2*(q(3)*q(4)-q(1)*q(2));
    C(3,1) = 2*(q(2)*q(4)-q(1)*q(3));
    C(3,2) = 2*(q(3)*q(4)+q(1)*q(2));
    C(3,3) = q(1)^2-q(2)^2-q(3)^2+q(4)^2;
end

function C = Quaternion2nCb(q)
    C = zeros(3);
    C(1,1) = q(1)^2+q(2)^2-q(3)^2-q(4)^2;
    C(2,1) = 2*(q(2)*q(3)-q(1)*q(4));
    C(3,1) = 2*(q(2)*q(4)+q(1)*q(3));
    C(1,2) = 2*(q(2)*q(3)+q(1)*q(4));
    C(2,2) = q(1)^2-q(2)^2+q(3)^2-q(4)^2;
    C(3,2) = 2*(q(3)*q(4)-q(1)*q(2));
    C(1,3) = 2*(q(2)*q(4)-q(1)*q(3));
    C(2,3) = 2*(q(3)*q(4)+q(1)*q(2));
    C(3,3) = q(1)^2-q(2)^2-q(3)^2+q(4)^2;
end

function q = Euler2Quaternion(pitch, roll, yaw)
    q=[cos(yaw/2)*cos(pitch/2)*cos(roll/2)-sin(yaw/2)*sin(pitch/2)*sin(roll/2);
       cos(yaw/2)*sin(pitch/2)*cos(roll/2)-sin(yaw/2)*cos(pitch/2)*sin(roll/2);
       cos(yaw/2)*cos(pitch/2)*sin(roll/2)+sin(yaw/2)*sin(pitch/2)*cos(roll/2);
       cos(yaw/2)*sin(pitch/2)*sin(roll/2)+sin(yaw/2)*cos(pitch/2)*cos(roll/2)];
end
