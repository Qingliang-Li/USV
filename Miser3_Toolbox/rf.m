function [rlx,rly,k,rb]=rf(lx,ly)
%输入：内四边形数据
%输出：外四边形顶点坐标，斜率，以及截距
for i=1:4
    k(i)=(ly(i+1)-ly(i))/(lx(i+1)-lx(i));%斜率
    b(i)=ly(i)-k(i)*lx(i);%截距
end
delta=8;%安全距离
rb=[b(1)-delta b(2)-delta b(3)+delta b(4)+delta];
rb=[rb rb(1)];
k=[k k(1)];
%计算坐标
for i=1:4
    rlx(i)=(rb(i+1)-rb(i))/(k(i)-k(i+1));
    rly(i)=k(i)*rlx(i)+rb(i);
end
rlx=[rlx rlx(1)];
rly=[rly rly(1)];
end