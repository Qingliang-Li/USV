function f = ocf(t,x,u,z,upar)

% Return the value of the right hand side of the state equations
% as a column vector.
%无人船参数
m1=141.85;m2=197.75;m3=15.6;
X=45.6;Y=29.54;Z=10.7;
%动态方程
f=z(1)*[x(4)*cos(x(3));
    x(4)*sin(x(3));
    u(1);
    u(2)];
end
% [x(4)*cos(x(3))-x(5)*sin(x(3));
%     x(4)*sin(x(3))+x(5)*cos(x(3));
%     x(6);
%     (u(1)+m2*x(5)*x(6)-X*x(4))/m1;
%     (-m1*x(4)*x(6)-Y*x(5))/m2;
%     (u(2)+(m1-m2)*x(4)*x(5)-Z*x(6))/m3];
