%  给定各个关节的驱动力矩，验证正向动力学的正确性
clc;clear;close all
global Qi J_type BB SS SE S0
global cc c0 ce
global m0 inertia0 m inertia mass
global d_time
global Qe
global Gravity
global num_q Ez
% Unit vector in z-axis
     Ez = [0;0;1];               %代表Z轴的向量
% Link Parameters
%%%   需要下面这些量   R0, A0, v0, w0, q, qd, F0, T0, Fe, Te, tau 
    BB = [ 0 1 2 3 4 5 6 ];
    SS = [ -1  1  0  0  0  0  0;       %%  关联矩阵
                 0 -1  1  0  0  0  0;
                 0  0 -1  1 0  0   0;
                 0  0  0 -1  1  0  0;
                 0  0  0  0 -1  1  0;
                 0  0  0  0  0 -1  1;
                 0  0  0  0  0  0 -1 ];
    S0 = [ 1 0 0 0 0 0 0 ];%  与基座相邻的为1，其余为0。
    SE = [ 0 0 0 0 0 0 1 ]; %  SE 为与外力矩相关的一个矩阵   如果为末端杆件 则SE(i)=0;  对于主程序中给出的为双臂机械臂，因此本文为【0 0 1 0 0 1】
    J_type = [ 'R' 'R' 'R' 'R' 'R' 'R' 'R' ];   %  关节类型
    c0 = [ 0   0 0 0 0 0 0;         %  杆件0 的质心到各个关节的i 的距离 ;在计算杆件i的质心位置时用到了c0
              0    0 0 0 0 0 0;
              3.5 0 0 0 0 0 0];
     m0 =15;  % 基座质量
     inertia0 =0.001*[ 2400  0  0;   % 基座转动惯量
              0 2400  0;
              0  0 2400 ];
     m = [ 20.76  20.76 40.45 31.15 20.76 20.76 27.96];  %  杆件质量
     mass = sum(m) + m0;    % 系统总质量
     inertia = [  0.645  0   0         0.645  0   0        0.284  0   0          0.137  0  0        0.645  0   0             0.645  0   0        1.311 0 0;   %  6个连杆的转动惯量
                   0 0.152   0          0 0.152   0        0 1.312 0.54             0 1.011 0           0  0.645   0            0  0.645   0        0 1.311 0;
                   0 0     0.645        0 0     0.645      0   0.54 1.308           0   0 0.013           0 0 0.152              0 0 0.152          0   0 0.220 ];
    ce = [ 0   0 0 0 0 0 0;         %  杆件i的质心到末端点的向量
              0    0 0 0 0 0 0;
              0    0 0 0 0 0 0.3525];
    Qe = [ 0   0 0 0 0 0 0;       %  末端点的方向
              0    0 0 0 0 0 0;
              0    0 0 0 0 0 0];

cc = zeros(3,7,7);
cc(:,1,1) = 0.001*[ 0 168.41 0 ]';       % 从连杆i的质心到关节i的向量 ,不知道是不是在连体系下的坐标？？？？
cc(:,2,2) = 0.001*[  0 -168.41 0]';       
cc(:,3,3) = 0.001*[ -968.5 0 -23.7]';
cc(:,4,4) = 0.001*[ -1162.32 0 0 ]';
cc(:,5,5) = 0.001*[ 0 0 -261.6 ]';
cc(:,6,6) =0.001* [ 0 0 -261.6 ]';
cc(:,7,7) =0.001* [ 0 0 -347.5 ]';

cc(:,1,2) =0.001* [ 0 -261.59 0]';
cc(:,2,3) =0.001* [ 0 261.59 0 ]';
cc(:,3,4) = 0.001*[ 1111.5 0 273.3]';
cc(:,4,5) =0.001* [ 917.68 0 0 ]';
cc(:,5,6) =0.001* [ 0 0 168.4]';
cc(:,6,7) =0.001* [ 0 0 168.4]';
% Initialize variables
q   = zeros(7,1);          % 机械臂各个关节角度
% % q=[0;90;0;0;0;-90;0]*pi/180;
qd  = zeros(7,1);
qdd = zeros(7,1);
v0  = [ 0 0 0 ]';   % 基座质心速度
w0  = [ 0 0 0 ]';   
vd0 = [ 0 0 0 ]';  %基座质心加速度 
wd0 = [ 0 0 0 ]';
vv = zeros(3,7);  % 各个杆件质心的速度
ww = zeros(3,7);
vd = zeros(3,7);
wd = zeros(3,7);
R0 = [ 0 0 0 ]';       %  基座的位置向量
Q0 = [ 0 0 0 ]';        % 基座RPY变换矩阵？？？？？？？？？
A0 = eye(3,3);       % 基座方向余弦
Fe = zeros(3,7);  %  作用在机械臂末端点上的力Fe和力矩Te
Te = zeros(3,7);
F0 = [ 0 0 0 ]';  %  作用在基座质心上的力F0和力矩T0
T0 = [ 0 0 0 ]';
tau = zeros(7,1);
d_time = 0.01;
Gravity = [ 0 0 0 ]';
num_q = length( q );
% 初始化 Qi 数组
Qi = zeros(3, num_q);
%%% Simulation Loop start
numj=1;
for time = 0:d_time:10,
        tau = zeros(7,1);
        tau(1:7,1)=[sin(time); 0.5*cos(time); 0.5*sin(time);0.1*cos(time); 0.1*sin(time); 0.05*sin(time); 0.01*sin(time)];
 AA = calc_aa( A0, q );  %%  各个杆件的方向  转换矩阵 
          RR = calc_pos( R0, A0, AA, q );   %%  杆件i的质心位置    
    [ R0, A0, v0, w0, q, qd ] = f_dyn_rk2( R0, A0, v0, w0, q, qd, F0, T0, Fe, Te, tau );
    v0j(numj,1:3)=v0;
    w0j(numj,1:3)=w0;
    numq(numj,1:7)=q;
     numj=numj+1;
end
 time = 0:d_time:10
 numqq=numq*180/pi;
    figure(111)
subplot(421)
    plot(time, numqq(:,1) ,'r-','LineWidth',1.8);grid on;legend('q1')
subplot(422)
    plot(time,numqq(:,2),'b-','LineWidth',1.8);grid on;legend('q2')
subplot(423)
    plot(time,numqq(:,3),'y-','LineWidth',1.8);grid on;legend('q3')
subplot(424)
    plot(time,numqq(:,4),'m-','LineWidth',1.8);grid on;legend('q4')
subplot(425)
    plot(time,numqq(:,5),'c-','LineWidth',1.8);grid on;legend('q5')
subplot(426)
    plot(time,numqq(:,6),'g-','LineWidth',1.8);grid on;legend('q6')
subplot(427)
    plot(time,numqq(:,7),'g-','LineWidth',1.8);grid on;legend('q7')
