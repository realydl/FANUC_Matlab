clc,clear;
%%
syms theta1 theta2 theta3 theta4 theta5 theta6;
syms dtheta1 dtheta2 dtheta3 dtheta4 dtheta5 dtheta6;
syms ddtheta1 ddtheta2 ddtheta3 ddtheta4 ddtheta5 ddtheta6;

qv=[theta1 theta2 theta3 theta4 theta5 theta6];
dqv=[dtheta1 dtheta2 dtheta3 dtheta4 dtheta5 dtheta6];
ddqv=[ddtheta1 ddtheta2 ddtheta3 ddtheta4 ddtheta5 ddtheta6];

l=[0.33 0.26 0.02 0.29 0.07];
m=[6 4 4 2 2 2];
r0=0.02;
%%
%���α���Ծ���
Ixx=[m(1)*l(1)^2/3 0 m(3)*l(4)^2/3 m(4)*l(5)^2/3 m(5)*l(5)^2/3 m(6)*l(5)^2/3];
Iyy=[m(1)*l(1)^2/3 m(2)*l(2)^2/3 m(3)*l(3)^2/3 m(4)*l(5)^2/3 0 m(6)*l(5)^2/3];
Izz=[0 m(2)*l(2)^2/3 m(3)*(l(3)^2+l(4)^2)/3 0 m(5)*l(5)^2/3 m(6)*r0^2/2];
Ixy=[0 0 -m(3)*l(3)*l(4)/3 0 0 0];
Ixz=zeros(1,6);
Iyz=zeros(1,6);
mx=[0 m(2)*l(2)/2 m(3)*l(3)/2 0 0 0];
my=[0 0 -m(3)*l(4)/2 0 -m(5)*l(5)/2 0];
mz=[-m(1)*l(1)/2 0 0 m(4)*l(5)/2 0 m(6)*l(5)/2];
for i=1:6
    J(:,:,i)=[(-Ixx(i)+Iyy(i)+Izz(i))/2 Ixy(i) Ixz(i) mx(i);
                Ixy(i) (Ixx(i)-Iyy(i)+Izz(i))/2 Iyz(i) my(i);
                Ixz(i) Iyz(i) (Ixx(i)+Iyy(i)-Izz(i))/2 mz(i);
                mx(i) my(i) mz(i) m(i)];
end
%%
%�������˸��ؽڱ任����T
T01=[cos(qv(1)) -sin(qv(1)) 0 0;
    sin(qv(1)) cos(qv(1)) 0 0;
    0 0 1 l(1)
    0 0 0 1];
T12=[cos(qv(2)) -sin(qv(2)) 0 0;
    0 0 -1 0;
    sin(qv(2)) cos(qv(2)) 0 0;
    0 0 0 1];
T23=[cos(qv(3)) -sin(qv(3)) 0 l(2);
    sin(qv(3)) cos(qv(3)) 0 0;
    0 0 1 0
    0 0 0 1];
T34=[cos(qv(4)) -sin(qv(4)) 0 l(3);
    0 0 -1 -l(4);
    sin(qv(4)) cos(qv(4)) 0 0
    0 0 0 1];
T45=[cos(qv(5)) -sin(qv(5)) 0 0;
    0 0 -1 0;
    sin(qv(4)) cos(qv(4)) 0 0;
    0 0 0 1];
T56=[cos(qv(6)) -sin(qv(6)) 0 0;
    0 0 -1 0;
    sin(qv(6)) cos(qv(6)) 0 0;
    0 0 0 1];

T01=T01;
T02=T01*T12;
T03=T01*T12*T23;
T04=T01*T12*T23*T34;
T05=T01*T12*T23*T34*T45;
T06=T01*T12*T23*T34*T45*T56;
T=cat(3,T01,T12,T23,T34,T45,T56);

%%
%����������D(q)
D=sym(zeros(6,6));
for i=1:6
    for j=1:6
        p=max(i,j);
        for n=p:6
            D(i,j)=D(i,j)+trace(diff(T(:,:,n),qv(j))...
                *J(:,:,n)*transpose(diff((T(:,:,n)),qv(i))));
        end
    end 
end
%%
%����������������ʸ��H
H=sym(zeros(6,1));
C=sym(zeros(6,6,6));
for i=1:6
    for k=1:6
        p=max([i j k]);
        for n=p:6
            C(j,k,i)=C(j,k,i)+trace(diff(diff(T(:,:,n),qv(j)),qv(k))...
            *J(:,:,n)*transpose(diff((T(:,:,n)),qv(i))));
        end
        H(i)=H(i)+C(j,k,i)*dqv(j)*dqv(k);
    end
end
%%
%�������ʸ��G
g=9.8;
gv=[0 0 g 0]';%���ٶ�����
%���������ڵ�ǰ��������ϵ�е�λ��
r=[0 0 -l(1)/2 1;
    l(2)/2 0 0 1;
    l(3)/2 -l(4)/2 0 1;
    0 0 l(5)/2 1;
    0 -l(5)/2 0 1;
    0 0 l(5)/2 1]';
G=sym(zeros(6,1));
for i=1:6
    for p=i:6
        G(i)=G(i)-m(p)*gv'*diff(T(:,:,p),qv(i))*r(:,p);
    end
end
%%
%��⶯��ѧ����
tauv=D*ddqv'+H+G;
%%
%FANUC �����˵Ĺ켣�滮
q0=deg2rad([0 90 0 0 0 0]);
qf=deg2rad([140 120 90 110 100 180]);
t=0:0.05:5;
[q,dq,ddq]=jtraj(q0,qf,t);%�õ������˸��ؽ�λ�á��ٶȺͼ��ٶ�
%%
%���FANUC�����˵Ĺؽ���������ʸ��M(t)
for i=1:length(t)
    %�ؽ�λ�ø�ֵ
    theta1=q(i,1);theta2=q(i,2);theta3=q(i,3);
    theta4=q(i,4);theta5=q(i,5);theta6=q(i,6);
    %�ؽ��ٶȸ�ֵ
    dtheta1=dq(i,1);dtheta2=dq(i,2);dtheta3=dq(i,3);
    dtheta4=dq(i,4);dtheta5=dq(i,5);dtheta6=dq(i,6);
    %�ؽڼ��ٶȸ�ֵ
    ddtheta1=ddq(i,1);ddtheta2=ddq(i,2);ddtheta3=ddq(i,3);
    ddtheta4=ddq(i,4);ddtheta5=ddq(i,5);ddtheta6=ddq(i,6);
    %����ʱ���ؽ���������
    tau(:,i)=eval(tauv);
end
%%
figure('Name','FANUC����������-�ؽڱ�����ϵ����');
%�ؽ�1
subplot(2,3,1);
plot(q(:,1),tau(1,:));
grid on
xlabel('theta_1(rad)');ylabel('�ؽ�1���أ�N*m��');
%�ؽ�2
subplot(2,3,2);
plot(q(:,2),tau(2,:));
grid on
xlabel('theta_2(rad)');ylabel('�ؽ�2���أ�N*m��');
%�ؽ�3
subplot(2,3,3);
plot(q(:,3),tau(3,:));
grid on
xlabel('theta_3(rad)');ylabel('�ؽ�3���أ�N*m��');
%�ؽ�4
subplot(2,3,4);
plot(q(:,4),tau(4,:));
grid on
xlabel('theta_4(rad)');ylabel('�ؽ�4���أ�N*m��');
%�ؽ�5
subplot(2,3,5);
plot(q(:,5),tau(5,:));
grid on
xlabel('theta_5(rad)');ylabel('�ؽ�5���أ�N*m��');
%�ؽ�6
subplot(2,3,6);
plot(q(:,6),tau(6,:));
grid on
xlabel('theta_6(rad)');ylabel('�ؽ�6���أ�N*m��');
    
    
    