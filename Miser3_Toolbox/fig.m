close all
clc
clear
set(0,'defaultfigurecolor','w')
%��Բ
theta=0:0.1:2*pi+0.1;
x0=[452 88 358 497 524;%x����
    436 120 78 306 155;%y
    34 26 34 45 53];%r��Բ�뾶
pho=x0(3,:)-8;%��Բ�뾶
%%���ı���
lx=[162,198,380,244,162];%��������
ly=[80,60,185,205,80];
[rlx,rly,k,rb]=rf(lx,ly);%�������ıߵ�����

lx1=[127,210,351,234,127];
ly1=[274,245,356,391,274];
[rlx1,rly1,k1,rb1]=rf(lx1,ly1);
%xx(1)=200;xx(2)=300;
%max(0,max(0,xx(2)-k1(1)*xx(1)-rb1(1))*min(0,xx(2)-k1(3)*xx(1)-rb1(3))*max(0,xx(2)-k1(2)*xx(1)-rb1(2))*min(0,xx(2)-k1(4)*xx(1)-rb1(4)))
figure(1)
p1=plot(0,0,'ro',500,500,'rp');%��ʼ��
hold on
for i=1:size(x0,2)
    %��Բ
    x(i,:)=x0(1,i)+x0(3,i)*cos(theta);
    y(i,:)=x0(2,i)+x0(3,i)*sin(theta);
    plot(x(i,:),y(i,:),'k--','HandleVisibility','off')
    %����Բ�����
    rx(i,:)=x0(1,i)+pho(i)*cos(theta);
    ry(i,:)=x0(2,i)+pho(i)*sin(theta);
    patch(transpose(rx(i,:)),transpose(ry(i,:)),'b','HandleVisibility','off')
    %fill(x(i,:),y(i,:));
end
%���ı���
plot(rlx,rly,'k--','HandleVisibility','off')
patch(lx,ly,'b','HandleVisibility','off')
plot(rlx1,rly1,'k--','HandleVisibility','off')
patch(lx1,ly1,'b','HandleVisibility','off')
load('x.mat');%���غ�������
p3=plot(xqd(1,:),xqd(2,:),'r','linewidth',1.2);%������
xlabel('{\itx}/m');
ylabel('{\ity}/m');
legend('���','�յ�','����')
grid on
set(gca,'GridLineStyle','--','GridColor','k','GridAlpha',0.1);
%legend(p2,'����')


n=size(xqd,2);
t=98.26;
dt=0:1/(n-1):1;
dt=dt*t;
Vmax=8;Vmin=0;
%���ٶ���ʱ���ͼ
figure(2)
plot(dt,xqd(4,:),'r','linewidth',1.2)
hold on
plot([0,dt(end)],[Vmax,Vmax],'k--',[0,dt(end)],[Vmin,Vmin],'k--')
axis([-inf 100,-1,9])
xlabel('{\itt}/s');
ylabel('{\fontname{Times New Roman}\it V}/(m��s^{-1})');
legend('{\fontname{Times New Roman}\it V}({\itt})','�߽�ֵ')

%���ؿ�����
load('u.mat')
amax=2;amin=-2;
rmax=0.1;rmin=-0.1;
%����������ʱ���ͼ
figure(3)
subplot(2,1,1)
plot(dt,uqd(1,:),'k','linewidth',1.2)
hold on
plot([0,dt(end)],[rmax,rmax],'k--',[0,dt(end)],[rmin,rmin],'k--')
axis([-inf 100,-0.15,0.15])
xlabel('{\itt}/s');
ylabel('���ٶ�/(rad��s^{-1})');
legend('{\it \gamma} ({\itt})','�߽�ֵ')
subplot(2,1,2)
plot(dt,uqd(2,:),'k','linewidth',1.2)
hold on
plot([0,dt(end)],[amax,amax],'k--',[0,dt(end)],[amin,amin],'k--')
axis([-inf 100,-2.5,2.5])
xlabel('{\itt}/s');
ylabel('���ٶ�/(m��{s^{-2}})');
legend('{\fontname{Times New Roman}\it{a}} ({\itt})','�߽�ֵ')
%plot(x,y,'r^')
%��xy��ʱ���ͼ
figure(4)
subplot(2,1,1)
plot(dt,xqd(1,:),'k','linewidth',1.2)
xlabel('{\itt}/s');
ylabel('{\itx}/m');
subplot(2,1,2)
plot(dt,xqd(2,:),'k','linewidth',1.2)
xlabel('{\itt}/s');
ylabel('{\ity}/m');

%�������ı������ݵ��Ӻ���
function [rlx,rly,k,rb]=rf(lx,ly)
for i=1:4
    k(i)=(ly(i+1)-ly(i))/(lx(i+1)-lx(i));
    b(i)=ly(i)-k(i)*lx(i);
end
delta=8;
rb=[b(1)-delta b(2)-delta b(3)+delta b(4)+delta];
rb=[rb rb(1)];
k=[k k(1)];
for i=1:4
    rlx(i)=(rb(i+1)-rb(i))/(k(i)-k(i+1));
    rly(i)=k(i)*rlx(i)+rb(i);
end
rlx=[rlx rlx(1)];
rly=[rly rly(1)];
end
