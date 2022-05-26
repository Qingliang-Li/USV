function phi = ocphi(ig,taut,x,z,upar,rhoab,kabs)
phi=0;
alpha=-1.5;
beta=3;
delta=1e8;
if ig == 0, phi=z(1)+delta*z(2)^beta+z(2)^alpha*((x(1)-500)^2+(x(2)-500)^2); end%时间T
if ig == 1, phi=x(1)-500; end%x坐标
if ig == 2, phi=x(2)-500; end%y坐标