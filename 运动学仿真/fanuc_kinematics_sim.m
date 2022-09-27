clc,clear;
%%
%建立FANUC机器人仿真模型
l=[0.33 0.26 0.02 0.29];
L(1)=Link('d',l(1),'a',0,'alpha',0,'modified','qlim',[-170 170]*pi/180);
L(2)=Link('d',0,'a',0,'alpha',pi/2,'modified','qlim',[-115 115]*pi/180);
L(3)=Link('d',0,'a',l(2),'alpha',0,'modified','qlim',[-201 201]*pi/180);
L(4)=Link('d',l(4),'a',l(3),'alpha',pi/2,'modified','qlim',[-190 190]*pi/180);
L(5)=Link('d',0,'a',0,'alpha',pi/2,'modified','qlim',[-180-120 180+120]*pi/180);
L(6)=Link('d',0,'a',0,'alpha',pi/2,'modified','qlim',[-360 360]*pi/180);
fanuc=SerialLink(L,'name','fanuc 200id');

qc=[0 pi/2 0 0 pi 0];
fanuc.plot(qc);
fanuc.teach(qc);
%%FANUC机器人工作空间
qL=fanuc.qlim;
N=90000;
ws_q=zeros(N,6);
a=rand(size(ws_q));

for j=1:N
    for i=1:6
        ws_q(j,i)=qL(i,1)+(qL(i,2)-qL(i,1))*a(j,i);
    end
end
%求机器人位置坐标
ws_x=cos(ws_q(:,1)).*(l(2)*cos(ws_q(:,2))+l(3)*cos(ws_q(:,2)+ws_q(:,3))+l(4)*sin(ws_q(:,2)+ws_q(:,3)));
ws_y=sin(ws_q(:,1)).*(l(2)*cos(ws_q(:,2))+l(3)*cos(ws_q(:,2)+ws_q(:,3))+l(4)*sin(ws_q(:,2)+ws_q(:,3)));
ws_z=l(1)+l(2)*sin(ws_q(:,2))+l(3)*sin(ws_q(:,2)+ws_q(:,3))-l(4)*cos(ws_q(:,2)+ws_q(:,3));

figure('Name','FANUC机器人工作空间')
plot3(ws_x,ws_y,ws_z,'r.');
grid;
xlabel('X');ylabel('Y');zlabel('Z');
%%
%初始点P0，逆运动学求解
P0=[0.3 0.25 0.6 pi/4 pi/2 -pi/4];
w0=P0(4);p0=P0(5);r0=P0(6);
R0=rotz(r0)*roty(p0)*rotx(w0);
T0=[R0 P0(1:3)';0 0 0 1];
%逆解theta0_1
theta0_1=atan2(T0(2,4),T0(1,4));
%逆解theta0_2
A0_2=((T0(1,4)*cos(theta0_1)+T0(2,4)*sin(theta0_1))^2+(T0(3,4)-l(1))^2+l(2)^2-l(3)^2-l(4)^2)...
/(2*l(2)*sqrt((T0(1,4)*cos(theta0_1)+T0(2,4)*sin(theta0_1))^2+(T0(3,4)-l(1))^2));
phi0_2=atan2((T0(1,4)*cos(theta0_1)+T0(2,4)*sin(theta0_1)),(T0(3,4)-l(1)));
theta0_2=atan2(A0_2,-sqrt(1-A0_2^2))-phi0_2;
%逆解theta0_3
A0_3=atan2(T0(1,4)*cos(theta0_1)+T0(2,4)*sin(theta0_1)-l(2)*cos(theta0_2),-T0(3,4)+l(1)+l(2)*sin(theta0_2));
phi0_3=atan2(l(3),l(4));
theta0_3=A0_3-theta0_2-phi0_3;
%逆解theta0_5
A0_5=T0(3,3)*cos(theta0_2+theta0_3)-(T0(1,3)*cos(theta0_1)+T0(2,3)*sin(theta0_1))*sin(theta0_2+theta0_3);
theta0_5=atan2(-sqrt(1-A0_5^2),A0_5);
%逆解theta0_4
theta0_4=atan2((T0(1,3)*sin(theta0_1)-T0(2,3)*cos(theta0_1))/sin(theta0_5),(T0(3,3)*sin(theta0_2+theta0_3)...
+(T0(1,3)*cos(theta0_1)+T0(2,3)*sin(theta0_1))*cos(theta0_2+theta0_3))/sin(theta0_5));
%逆解theta0_6
theta0_6=atan2((T0(3,2)*cos(theta0_2+theta0_3)-(T0(1,2)*cos(theta0_1)+T0(2,2)*sin(theta0_1))*sin(theta0_2+theta0_3))/sin(theta0_5),...
(T0(3,1)*cos(theta0_2+theta0_3)-(T0(1,1)*cos(theta0_1)+T0(2,1)*sin(theta0_1))*sin(theta0_2+theta0_3))/-sin(theta0_5));

q0=[theta0_1 theta0_2 theta0_3 theta0_4 theta0_5 theta0_6];

%终止点Pf，逆运动学求解
Pf=[0.1 0.4 0.4 0 0 0];
wf=Pf(4);
pf=Pf(5);
rf=Pf(6);
Rf=rotz(rf)*roty(pf)*rotx(wf);
Tf=[Rf Pf(1:3)';0 0 0 1];
%逆解theta0_1
thetaf_1=atan2(Tf(2,4),Tf(1,4));
%逆解theta0_2
Af_2=((Tf(1,4)*cos(thetaf_1)+Tf(2,4)*sin(thetaf_1))^2+(Tf(3,4)-l(1))^2+l(2)^2-l(3)^2-l(4)^2)...
/(2*l(2)*sqrt((Tf(1,4)*cos(thetaf_1)+Tf(2,4)*sin(thetaf_1))^2+(Tf(3,4)-l(1))^2));
phif_2=atan2((Tf(1,4)*cos(thetaf_1)+Tf(2,4)*sin(thetaf_1)),(Tf(3,4)-l(1)));
thetaf_2=atan2(Af_2,-sqrt(1-Af_2^2))-phif_2;
%逆解theta0_3
Af_3=atan2(Tf(1,4)*cos(thetaf_1)+Tf(2,4)*sin(thetaf_1)-l(2)*cos(thetaf_2),-Tf(3,4)+l(1)+l(2)*sin(thetaf_2));
phif_3=atan2(l(3),l(4));
thetaf_3=Af_3-thetaf_2-phif_3;
%逆解theta0_5
Af_5=Tf(3,3)*cos(thetaf_2+thetaf_3)-(Tf(1,3)*cos(thetaf_1)+Tf(2,3)*sin(thetaf_1))*sin(thetaf_2+thetaf_3);
thetaf_5=atan2(-sqrt(1-Af_5^2),Af_5);
%逆解theta0_4
thetaf_4=atan2((Tf(1,3)*sin(thetaf_1)-Tf(2,3)*cos(thetaf_1))/sin(thetaf_5),(Tf(3,3)*sin(thetaf_2+thetaf_3)...
+(Tf(1,3)*cos(thetaf_1)+Tf(2,3)*sin(thetaf_1))*cos(thetaf_2+thetaf_3))/sin(thetaf_5));
%逆解theta0_6
thetaf_6=atan2((Tf(3,2)*cos(thetaf_2+thetaf_3)-(Tf(1,2)*cos(thetaf_1)+Tf(2,2)*sin(thetaf_1))*sin(thetaf_2+thetaf_3))/sin(thetaf_5),...
(Tf(3,1)*cos(thetaf_2+thetaf_3)-(Tf(1,1)*cos(thetaf_1)+Tf(2,1)*sin(thetaf_1))*sin(thetaf_2+thetaf_3))/-sin(thetaf_5));

qf=[thetaf_1 thetaf_2 thetaf_3 thetaf_4 thetaf_5 thetaf_6];

%%
%规划轨迹
t=0:0.1:5;
[q qd qdd]=jtraj(q0,qf,t);
%正运动学仿真
T=double(fanuc.fkine(q));
[x y z]=transl(T);

%绘制曲线图
figure('Name','“FANUC Robot”机器人正运动学演示');
fanuc.plot(q);
figure('Name','FANUC 机器人末端运动轨迹');
plot3(x,y,z,'r-o','MarkerFaceColor','r');
grid;
xlabel('X');ylabel('Y');zlabel('Z');
figure('Name','各关节的位移速度加速度曲线');
%关节1
subplot(3,6,1);
plot(t,q(:,1));
grid on;
xlim([0 t(end)]);
title('关节1');
ylabel('位移（rad）');
subplot(3,6,7);
plot(t,qd(:,1));
grid on;
xlim([0 t(end)]);
ylabel('速度（rad/s）');
subplot(3,6,13);
plot(t,qdd(:,1));
grid on;
xlim([0 t(end)]);
xlabel('时间（s）');ylabel('加速度（rad/^2s）');
%关节2
subplot(3,6,2);
plot(t,q(:,2));
grid on;
xlim([0 t(end)]);
title('关节2');
subplot(3,6,8);
plot(t,qd(:,2));
grid on;
xlim([0 t(end)]);
subplot(3,6,14);
plot(t,qdd(:,2));
grid on;
xlim([0 t(end)]);
xlabel('时间（s）');
%关节3
subplot(3,6,3);
plot(t,q(:,3));
grid on;
xlim([0 t(end)]);
title('关节3');
subplot(3,6,9);
plot(t,qd(:,3));
grid on;
xlim([0 t(end)]);
subplot(3,6,15);
plot(t,qdd(:,3));
grid on;
xlim([0 t(end)]);
xlabel('时间（s）');
%关节4
subplot(3,6,4);
plot(t,q(:,4));
grid on;
xlim([0 t(end)]);
title('关节4');
subplot(3,6,10);
plot(t,qd(:,4));
grid on;
xlim([0 t(end)]);
subplot(3,6,16);
plot(t,qdd(:,4));
grid on;
xlim([0 t(end)]);
xlabel('时间（s）');
%关节5
subplot(3,6,5);
plot(t,q(:,5));
grid on;
xlim([0 t(end)]);
title('关节5');
subplot(3,6,11);
plot(t,qd(:,5));
grid on;
xlim([0 t(end)]);
subplot(3,6,17);
plot(t,qdd(:,5));
grid on;
xlim([0 t(end)]);
xlabel('时间（s）');
%关节6
subplot(3,6,6);
plot(t,q(:,6));
grid on;
xlim([0 t(end)]);
title('关节6');
subplot(3,6,12);
plot(t,qd(:,6));
grid on;
xlim([0 t(end)]);
subplot(3,6,18);
plot(t,qdd(:,6));
grid on;
xlim([0 t(end)]);
xlabel('时间（s）');


