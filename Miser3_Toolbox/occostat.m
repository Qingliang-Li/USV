function dpsi = occostat(t,psi,flag,zt,kpright,kpleft,kint,kpar,...
 koncts,tk,kptk,kpcp,tqd,lcint,lget,ngc,epsjt,rhoab,kabs,upar,numder,reltest)

global Cb
global xqd dxqd

% obtain the value of u
ut=[];
for ic=1:length(koncts)
 icp=kpar(ic,kint);
 if koncts(ic)==0
  ut(ic)=Cb(icp);
 else
  leftindex=kptk(ic)+icp-kpcp(ic);
  tkleft=tk(leftindex);
  tkright=tk(leftindex+1);
  ut(ic)=Cb(icp)+(Cb(icp+1)-Cb(icp))*(t-tkleft)/(tkright-tkleft);
 end
end
ut=ut';
	
% given t, determine 3 points for the Hermite interpolation on x
% kherm is the leftmost of the 3 indices
ind=find(tqd>=t);
if t==tqd(ind(1))
 xt=xqd(:,ind(1));
else
 kherm=ind(1)-1;
 if kherm < 1, kherm=1; end
 if kherm == kpright, kherm=kherm-1; end
 if (abs(t-tqd(kherm)) < abs(t-tqd(kherm+1))), kherm=kherm-1; end
 if kherm == kpright-1, kherm=kherm-1; end
 if kherm == kpleft-1, kherm=kpleft; end

% Do Hermite interpolation
 tk=tqd(kherm:kherm+2);
 h=tk(2)-tk(1);
 tt=t*ones(1,3)-tk;
 tt2=tt.*tt;
 l2=[0.25*tt2(2)*tt2(3) tt2(1)*tt2(3) 0.25*tt2(1)*tt2(2)];
 w=[(1+3*tt(1)/h) 1 (1-3*tt(3)/h)].*l2;
 wd=tt.*l2;
 xt =xqd(:,kherm:kherm+2)*w';
 if kherm+2 < kpright
  xt=xt+dxqd(:,kherm:kherm+2)*wd';
 else
%  make sure dxdt is the limit from the left of kpright.
  df = ocf(tqd(kpright),xqd(:,kherm+2),ut,zt,upar);
  xt=xt+[dxqd(:,kherm:kherm+1) df]*wd';
 end
 xt = xt./h^4;
end

ns=length(xt);
if ~numder(1,2)
 dg0 = ocdg0dx(t,xt,ut,zt,upar,rhoab,kabs);
else
 dg0 = cnumdg0(t,xt,ut,zt,upar,rhoab,kabs,reltest);
end
if ~numder(1,1)
 df = ocdfdx(t,xt,ut,zt,upar);
else
 df = cnumdf(t,xt,ut,zt,upar,reltest);
end
dpsi = -dg0'-df'*psi(1:ns,1);

if ngc > 0
 g = ocg(t,xt,ut,zt,upar,rhoab,kabs);
 if ~numder(1,4)
  dg = ocdgdx(t,xt,ut,zt,upar,rhoab,kabs);
 else
  dg = cnumdg(t,xt,ut,zt,upar,rhoab,kabs,reltest);
 end
 for ig=1:ngc
  if lcint(ig)
   if lget(ig)
    g(ig)=ocdetsm(g(ig),epsjt);
   else
    g(ig)=1;
   end
   dpdt = -dg(ig,:)*g(ig)-psi(ig*ns+1:ns*(ig+1),1)'*df;
   dpsi = [dpsi; dpdt'];
  else
   dpsi = [dpsi zeros(ns,1)];
  end
 end
end
% End of occostat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dgeps = ocdetsm(g,epsjt)
% Derivative wrt argument of ocetsm.

if g >= epsjt
 dgeps=0;
elseif g <= -epsjt
 dgeps=1;
else
 dgeps=-0.5*(g-epsjt)/epsjt;
end

% End of ocdetsm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dfdx = cnumdf(t,x,u,z,upar,reltest);
% Finite difference code for derivative of f wrt x. Uses central differences.
% Not recommended for very accurate work, and it is slow.

ha=eps^(1/3);
nx=length(x);
for i=1:nx
	xi=x(i);
	hh=ha*max(reltest,abs(xi));
	x(i)=xi+hh;
	fph = ocf(t,x,u,z,upar);
	x(i)=xi-hh;
	fmh = ocf(t,x,u,z,upar);
	dfdx(:,i)=(fph-fmh)/(2*hh);
	x(i)=xi;
end

% End of cnumdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&
function dg0dx = cnumdg0(t,x,u,z,upar,rhoab,kabs,reltest);
% Finite difference code for derivative of g0 wrt x. Uses central differences.
% Not recommended for very accurate work, and it is slow.

ha=eps^(1/3);
nx=length(x);
for i=1:nx
	xi=x(i);
	hh=ha*max(reltest,abs(xi));
	x(i)=xi+hh;
	g0ph = ocg0(t,x,u,z,upar,rhoab,kabs);
	x(i)=xi-hh;
	g0mh = ocg0(t,x,u,z,upar,rhoab,kabs);
	dg0dx(i)=(g0ph-g0mh)/(2*hh);
	x(i)=xi;
end

% End of cnumdg0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&
function dgdx = cnumdg(t,x,u,z,upar,rhoab,kabs,reltest);
% Finite difference code for derivatives of g wrt x. Uses central differences.
% Not recommended for very accurate work, and it is slow.

ha=eps^(1/3);
nx=length(x);
for i=1:nx
	xi=x(i);
	hh=ha*max(reltest,abs(xi));
	x(i)=xi+hh;
	gph = ocg(t,x,u,z,upar,rhoab,kabs);
	x(i)=xi-hh;
	gmh = ocg(t,x,u,z,upar,rhoab,kabs);
	dgdx(:,i)=(gph-gmh)/(2*hh);
	x(i)=xi;
end

% End of cnumdg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&

