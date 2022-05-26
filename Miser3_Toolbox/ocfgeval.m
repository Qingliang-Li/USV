function [f0,gradf0,g,gradg]=ocfgeval(iflag,c,tolx,reltest,kpqk,tknot,...
 dtk,tk,tqd,tau,kpqtau,kpar,koncts,kpcp,kptk,ns,nz,ngc,lget,numder,...
 epsjt,taujt,rhoab,tolpsi,tfinal,ferr,kout,misc)

%  Formal parameters.
%  iflag indicates which function values to return.
%        =0001  f0,  =0010  g,  =0100  gradf0,  =1000  gradg
%   if any digit is 2, compute even if same c, and also recompute
%   relevant differential equations if necessary.

%  array cspar,  value of control parameters and system parameters.
%  scalar f0, value of canonical objective function.
%  array g,  value of the canonical constraints.
%  array gradf0,  value of gradient of f0 with respect to parameters.
%  array gradg,  value of gradients of constraints.

%  The name of the objective function has been changed to g0 for the
%  user but in this set of subroutines we use f0 so that names of
%  variables are different in the first six characters.

%  Calculates the canonical optimal control problem functionals and
%  gradients. This routine is meant to be the internal engine of the
%  package as all function values and gradient values for the 
%  optimization routine will be evaluated here. It is envisaged that the
%  routines required by FMINCON will call this routine with
%  iflag set appropriately. This routine will keep the values of the
%  last control parameters on exit so that on re-entry the values of the
%  control parameters can be checked to see if they have changed. In
%  this way duplication of evaluation can be avoided. Different
%  optimization routines differ as to whether objective function or
%  constraints are evaluated first and whether gradients are computed 
%  separately or not.
%  This routine calls a MATLAB numerical integration of ode routine.

%  The numerical integration and quadrature are computed over each
%  subinterval defined by all knot sets and internal points of the
%  problem. 

%_____________________________________________________________________________
%  logical variables

%  lsame=.true. indicates same values for parameters.
%  lcstate=.true. indicates state equations already integrated.
%  lccostate=.true. indicates costate equations already integrated.
%  lcf0=.true. indicates objective function already computed.
%  lcg=.true. indicates constraints already computed.
%  lcgradf0=.true. indicates gradient of objective already computed.
%  lcgradg=.true. indicates gradient of constraints already computed.
%  lstate=.true. indicates state equations need to be integrated.
%  lcostate=.true. indicates costate equations need to be integrated.
%  lf0=.true. indicates objective function needs to be computed.
%  lg=.true. indicates constraints needs to be computed.
%  lgradf0=.true. indicates gradient of objective needs to be computed.
%  lgradg=.true. indicates gradient of constraints needs to be computed.

global Cb Logicvalues Ntimes
global xqd dxqd uqd psiqd
% Cb is initially set to 999*ones(1,nctot)
% Logicvalues=[lsame,lstate,lccostate,lcf0,lcg,lcgradf0,lcgradg]
% ntimes is the number of times ocfgeval is called
% these are initially set to zero in the first call to ocfgeval from the main
%_____________________________________________________________________________

f0=[];gradf0=[];g=[];gradg=[];

% Extract some parameters
nctot=length(c);
nknot=length(tknot);
nps=ns*(ngc+1);
nct=nctot-nz;
nqp=length(tqd);
if isempty(koncts)
 nc=0;
else
 nc=length(koncts);
end
tstart=tqd(1);
ncheck=misc(1);
kabs=misc(5);
upar=misc(7:6+misc(6));

lsame=Logicvalues(1);
lstate=Logicvalues(2);
lccostate=Logicvalues(3);
lcf0=Logicvalues(4);
lcg=Logicvalues(5);
lcgradf0=Logicvalues(6);
lcgradg=Logicvalues(7);

% Check same set of control parameters and variables.
zt=Cb(nct+1:end);
lsame=isequal(c,Cb);
if ~lsame | iflag==2222	 
% New values for controls and variables.
 lcstate=0;
 lccostate=0;
 lcf0=0;
 lcg=0;
 lcgradf0=0;
 lcgradg=0;
 lsame=0;
 zt=c(nct+1:end);
 Cb=c;
end

% Sort out what is required from what is already computed.
fl4=rem(iflag,10);
fl3=rem((iflag-fl4)/10,10);
fl2=rem((iflag-10*fl3-fl4)/100,10);
fl1=rem((iflag-100*fl2-10*fl3-fl4)/1000,10);
if (fl4==2 | fl3==2) lcstate=0; end
lf0=((fl4==1 & ~lcf0) | fl4==2);
lg=((fl3==1 & ~lcg) | fl3==2);
lstate=(lf0 | lg) & ~lcstate;
lgradf0=((fl2==1 & ~lcgradf0) | fl2==2);
lgradg=((fl1==1 & ~lcgradg) | fl1==2);
if (fl2==2 | fl1==1), lccostate=0; end
lcostate=((lgradf0 | lgradg) & ~lccostate);
if ngc==0
 lg=0;
 lgradg=0;
end

%_____________________________________________________________________________		
if lstate
% Integrate the state equations between knots
 options=odeset('RelTol',tolx,'AbsTol',tolx*reltest,'Jacobian','off');
 xt = ocxzero(zt,upar); % put initial values of x in xt
	
 for ik=1:nknot-1
  kint=ik;	% kint is a pointer to the current knot subinterval
  kpleft=kpqk(ik);
  kpright=kpqk(ik+1);
  tsteps=tqd(kpleft:kpright);
% ODE solver
  if misc(2)==0
   [tde,yde,stats]=ode45('ocstate',tsteps,xt,options,zt,kint,...
    kpar,koncts,tk,kptk,kpcp,upar);
	if stats(3)>1000 & length(yde) == kpright-kpleft+1
    fprintf('\n state odes are probably stiff, you may wish to restart');
    fprintf('with the sstate stiff option set to ''1''\n');
   end
  else
   [tde,yde,stats]=ode15s('ocstate',tsteps,xt,options,zt,kint,...
    kpar,koncts,tk,kptk,kpcp,upar);
  end
 
% error in ode solver - write to error file
  lenyde=size(yde,1);
  if lenyde < kpright-kpleft+1
   iqq=kpleft+lenyde-1;	
   fprintf('\n\n State equations failed to integrate\n')
   fprintf(' Look to the error file for the state values up to this');
   fprintf(' time =%10.5f \n',tqd(iqq));
   fprintf(' The last state values computed are:\n');
   fprintf(' %12.6e %12.6e %12.6e %12.6e %12.6e',yde(lenyde,:));
   ut=[];
   if nc > 0
    ut=occalu(tqd(iqq),kint,kpar,koncts,tk,kptk,kpcp);
    fprintf('\n The control values at this time are:\n');
    fprintf('%12.6e %12.6e %12.6e %12.6e %12.6e',ut);
   end
   if nz > 0
    fprintf('\n The system parameters are:\n');
	 fprintf('%12.6e %12.6e %12.6e %12.6e %12.6e',zt);
   end	  
   fprintf(ferr,'\n State equations failed to integrate\n');
   fprintf(ferr,' A history of the state equation integration follows for');
   fprintf(ferr,' each quadrature point.\n The line of stars indicates');
   fprintf(ferr,' the end of this attempt at integration.\n Data after ');
   fprintf(ferr,'the stars are from the previous successful integration.\n');
   if nz > 0
    fprintf(ferr,' The system parameters are:\n');
 	 fprintf(ferr,'%12.6e %12.6e %12.6e %12.6e %12.6e',zt);
   end
   for ip=1:kpleft-1
    fprintf(ferr,'\n quadrature point no. %4i time =%10.5f',ip,tqd(ip));
    ut=uqd(:,ip);
    fprintf(ferr,'\n controls  %12.6e %12.6e %12.6e %12.6e %12.6e',ut);
    fprintf(ferr,'\n states  %12.6e %12.6e %12.6e %12.6e %12.6e',xqd(:,ip));
	 dxqd(:,ip)=ocf(tqd(ip),xt,ut,zt,upar);
    fprintf(ferr,'\n dx/dt %12.6e %12.6e %12.6e %12.6e %12.6e',dxqd(:,ip));
 end
    for ip=kpleft:iqq-1
    fprintf(ferr,'\n quadrature point no. %4i time =%10.5f',ip,tqd(ip));
    ut=[];
    if nc > 0
     ut=occalu(tqd(ip),kint,kpar,koncts,tk,kptk,kpcp);	
     fprintf(ferr,'\n controls  %12.6e %12.6e %12.6e %12.6e %12.6e',ut);
    end
    fprintf(ferr,'\n states  %12.6e %12.6e %12.6e %12.6e %12.6e',xqd(:,ip));
	 dxqd(:,ip)=ocf(tqd(ip),xt,ut,zt,upar);
    fprintf(ferr,'\n dx/dt %12.6e %12.6e %12.6e %12.6e %12.6e',dxqd(:,ip));
   end
   fprintf(ferr,'\n*********************************************');
   for ip=iqq+1:nqp
    fprintf(ferr,'\n quadrature point no. %4i time =%10.5f',ip,tqd(ip));
    fprintf(ferr,'\n controls  %12.6e %12.6e %12.6e %12.6e %12.6e',uqd(:,ip));
    fprintf(ferr,'\n states  %12.6e %12.6e %12.6e %12.6e %12.6e',xqd(:,ip));
    fprintf(ferr,'\n dx/dt %12.6e %12.6e %12.6e %12.6e %12.6e',dxqd(:,ip));
   end
% stop execution
   fclose('all');
   error('Stopping Execution of MISER - all open files have been closed');
		 
  else		 
   xqd(:,kpleft:kpright)=yde';		
   for iq=kpleft:kpright
    xt=xqd(:,iq);
    tt=tqd(iq);
    ut=[];
    if nc > 0
     ut=occalu(tt,kint,kpar,koncts,tk,kptk,kpcp);
     uqd(:,iq)=ut;
    end
    dxqd(:,iq) = ocf (tt,xt,ut,zt,upar);
   end
  end
 end
 lcstate=1;
end 

%_____________________________________________________________________________
% Check derivatives numerically if required.
reply='n';
analytd=0;
if any(~numder(1,1:4)), analytd=1; end
if any(~numder(2,1:2)) | numder(2,4)==0, analytd=1; end
if any(~numder(3,:)), analytd=1; end
if analytd & ncheck > 0 & rem(Ntimes,ncheck)==0
 worst=zeros(5,3);
 kflag=zeros(5,3);
 oworst=0;
 okflag=0;
 kint=1;
 indf0=zeros(1,ns+nc+nz);
 indf=zeros(ns,ns+nc+nz);
 indg=[];
 if ngc > 0, indg=zeros(ngc,ns+nc+nz); end
 indp=zeros(ngc+1,ns+nz);
 indx0=[];
 if nz > 0, indx0=zeros(ns,nz); end
 for iq=1:nqp
  if iq>=kpqk(kint+1), kint=min([kint+1,nknot-1]); end
  tt=tqd(iq);
  ut=[];
  if nc > 0, ut = occalu(tt,kint,kpar,koncts,tk,kptk,kpcp); end
  [worst,kflag,oworst,okflag,indf0,indf,indg,indp,indx0]=ocdcheck(tt,...
   xqd(:,iq),ut,zt,worst,kflag,oworst,okflag,indf0,indf,indg,indp,...
   indx0,ngc,tstart,tfinal,tau,rhoab,kabs,reltest,ferr,kout,upar,numder);
 end
 fprintf('\n\n Finished checking user derivatives.\n\n');
 if (kout(1)>=2 | okflag==2 | kout(1)*okflag==1)
  fprintf('##############################################');
  fprintf('##############################################\n');
  fprintf('Check on user supplied gradients against centred differences.');
  fprintf('\n okflag= %2i (0 = OK, 1 = fuzzy, 2 = wrong)\n',okflag);
  fprintf(' worst (actual - numerical) is %15.5e\n',oworst);
  fprintf('##############################################');
  fprintf('##############################################\n');
  fprintf(ferr,'\n\n ####################################################');
  fprintf(ferr,'##########################################\n');
  fprintf(ferr,'Check on user supplied gradients against centred differences.');
  fprintf(ferr,'\n okflag= %2i (0 = OK, 1 = fuzzy, 2 = wrong)\n',okflag);
  fprintf(ferr,' worst (actual - numerical) is %15.5e\n',oworst);
  fprintf(ferr,'####################################################');
  fprintf(ferr,'##########################################\n');
 end
 if okflag==2
  fprintf('\n Table of infinity norm of "actual - numerical" over all');
  fprintf(' quadrature points.\n The ok-fuzzy-wrong flag is also shown.');
  fprintf('\n\n Function   del /del x      del /del u      del /del z\n');
  fprintf('   f0     %11.4e%3i  %11.4e%3i  %11.4e%3i \n',...
   [worst(1,:); kflag(1,:)]);
  fprintf('   f      %11.4e%3i  %11.4e%3i  %11.4e%3i \n',...
   [worst(2,:); kflag(2,:)]);
  fprintf('   g      %11.4e%3i  %11.4e%3i  %11.4e%3i \n',...
   [worst(3,:); kflag(3,:)]);
  fprintf('   phi    %11.4e%3i                  %11.4e%3i \n',...
   worst(4,1),kflag(4,1),worst(4,3),kflag(4,3));
  fprintf('   x0                                     %11.4e%3i \n',...
   worst(5,3),kflag(5,3));
  fprintf('\n ##############################################');
  fprintf('##############################################\n');
  fprintf(ferr,'\n Table of infinity norm of "actual - numerical" over all');
  fprintf(ferr,' quadrature points.\n The ok-fuzzy-wrong flag is also shown.');
  fprintf(ferr,'\n\n Function   del /del x      del /del u      del /del z');
  fprintf(ferr,'\n   f0     %11.4e%3i  %11.4e%3i  %11.4e%3i \n',...
   [worst(1,:); kflag(1,:)]);
  fprintf(ferr,'   f      %11.4e%3i  %11.4e%3i  %11.4e%3i \n',...
   [worst(2,:); kflag(2,:)]);
  fprintf(ferr,'   g      %11.4e%3i  %11.4e%3i  %11.4e%3i \n',...
   [worst(3,:); kflag(3,:)]);
  fprintf(ferr,'   phi    %11.4e%3i                  %11.4e%3i \n',...
   worst(4,1),kflag(4,1),worst(4,3),kflag(4,3));
  fprintf(ferr,'   x0                                     %11.4e%3i \n',...
   worst(5,3),kflag(5,3));
  fprintf('##############################################');
  fprintf('##############################################\n');
  reply=input('\n Do you wish to abort this run? (y/n) ','s');
  if reply == 'y' 
   fprintf('\n Aborting this run due to wrong user derivatives as');
   fprintf('displayed above.\n Pattern of errors follow.\n');
   ociprint(ns,nc,nz,nqp,indf0,indf,indg,indp,indx0,kflag,ferr);
   fprintf('\n\n');
  end
 end
end
if reply == 'y'
 fclose('all');
 error('Stopping Execution of MISER - all open files have been closed');
end

%_____________________________________________________________________________
% Can now compute all function values at quadrature points.
% so compute the objective function and/or the constraints.
if lf0 | lg
% Compute the phi functions if necessary.
 if lf0, fj = ocphi(0,tfinal,xqd(:,nqp),zt,upar,rhoab,kabs); end	 
 if lg
  for ig=1:ngc
   fg(ig) = ocphi(ig,tau(ig),xqd(:,kpqtau(ig)),zt,upar,rhoab,kabs);
   if lget(ig),  fg(ig)=fg(ig)+taujt; end
  end
 end
% Compute the integrals - one sub-interval at a time.
 f0error=0;
 if ngc > 0, gerror=zeros(1,ngc); end
 for ik=1:nknot-1
  kint=ik;
  dt=dtk(ik);
  numints=kpqk(ik+1)-kpqk(ik);
% Compute Romberg weights
  w=[3.5 16 6 16];
  werror=[-0.25 1 -1.5 1];
  for i=2:floor(numints/4);
   w=[w 7 16 6 16];
   werror=[werror -0.5 1 -1.5 1];
  end
  w=4*dt*[w 3.5]/45;
  werror=4*dt*[werror -0.25]/45;
% Now compute integrands for objective and constraints
  fjint=[];gint=[];
  for iqq=1:numints+1
   iq=kpqk(ik)+iqq-1;
   tt=tqd(iq);
	xt=xqd(:,iq);
   ut=[];
   if nc > 0, ut = occalu(tt,kint,kpar,koncts,tk,kptk,kpcp); end
   if lf0
    fjint(iqq) = ocg0(tt,xt,ut,zt,upar,rhoab,kabs);
   end
   if lg
    gint(:,iqq) = ocg(tt,xt,ut,zt,upar,rhoab,kabs);
    for ig=1:ngc
     if lget(ig), gint(ig,iqq)=ocetsm(gint(ig,iqq),epsjt); end
     if kpqk(ik+1)>kpqtau(ig), gint(ig,iqq)=0; end
    end
   end
  end
% Now for Romberg computations and approx error.
  if lf0
   fj=fj+w*fjint';
   f0error=f0error+werror*fjint';
  end
  if lg
   fg=fg+w*gint';
   gerror=gerror+werror*gint';
  end
 end
 lcf0=max([lcf0,lf0]);
 lcg=max([lcg,lg]);
end

%_____________________________________________________________________________
% Done with function values - now for gradients.
if lcostate
% Integrate the costate equations, if necessary, over each sub-interval.
 options=odeset('RelTol',tolpsi,'AbsTol',tolpsi*reltest,'Jacobian','off');
 psit = zeros(nps,1); % set terminal values of costates
 for ik=nknot:-1:2
  kint=ik-1;	% kint is a pointer to the current knot subinterval
  kpleft=kpqk(kint);
  kpright=kpqk(ik);
  tsteps=tqd(kpright:-1:kpleft);
  [psit,lcint] = ocpsiset(xqd(:,kpright),zt,kpqk(nknot),ngc,kpright,...
   kpqtau,tau,tfinal,psit,rhoab,kabs,upar,numder,reltest);
 
% ODE solver
  if misc(3)==0
   [tde,yde,stats]=ode45('occostat',tsteps,psit,options,zt,...
    kpright,kpleft,kint,kpar,koncts,tk,kptk,kpcp,tqd,...
    lcint,lget,ngc,epsjt,rhoab,kabs,upar,numder,reltest);
   if stats(3)>1000 & length(yde) == kpright-kpleft+1
    fprintf('\n costate odes are probably stiff, you may wish to restart');
    fprintf('with the costate stiff option set to ''1''\n');
   end
  else
   [tde,yde,stats]=ode15s('occostat',tsteps,psit,options,zt,...
    kpright,kpleft,kint,kpar,koncts,tk,kptk,kpcp,tqd,...
    lcint,lget,ngc,epsjt,rhoab,kabs,upar,numder,reltest);				
  end	


 % error in ode solver - write to error file
  lenyde=size(yde,1);
  if lenyde < kpright-kpleft+1
   iqq=kpleft+lenyde-1;
   fprintf('\n\n Costate equations failed to integrate\n')
   fprintf(' Look to the error file for the state values and costate values');
   fprintf(' time =%10.5f \n',tqd(iqq));
   fprintf(' The last costate values computed are:');
	for ig=1:ngc+1
    fprintf('\n %12.6e %12.6e %12.6e %12.6e %12.6e',...
     yde(lenyde,(ig-1)*ns+1:ig*ns));
   end
   fprintf('\n The last state values used are:\n');
   fprintf(' %12.6e %12.6e %12.6e %12.6e %12.6e',xqd(:,iqq));
   ut=[];
   if nc > 0
    ut=occalu(tqd(iqq),kint,kpar,koncts,tk,kptk,kpcp);
    fprintf('\n The control values at this time are:\n');
    fprintf('%12.6e %12.6e %12.6e %12.6e %12.6e',ut);
   end
   if nz > 0
    fprintf('\n The system parameters are:\n');
	 fprintf('%12.6e %12.6e %12.6e %12.6e %12.6e',zt);
   end	  
   fprintf(ferr,'\n Costate equations failed to integrate\n');
   fprintf(ferr,' A history of the costate equation integration follows for');
   fprintf(ferr,' each quadrature point.\n The line of stars indicates');
   fprintf(ferr,' the end of this attempt at integration.\n Data after ');
   fprintf(ferr,'the stars are from the previous successful integration.\n');
   if nz > 0
    fprintf(ferr,' The system parameters are:\n');
 	 fprintf(ferr,'%12.6e %12.6e %12.6e %12.6e %12.6e',zt);
   end
   for ip=nqp:-1:iqq+1
    fprintf(ferr,'\n quadrature point no. %4i time =%10.5f',ip,tqd(ip));
    ut=[];
    if nc > 0
     ut=occalu(tqd(ip),kint,kpar,koncts,tk,kptk,kpcp);	
     fprintf(ferr,'\n controls  %12.6e %12.6e %12.6e %12.6e %12.6e',ut);
    end
    fprintf(ferr,'\n states  %12.6e %12.6e %12.6e %12.6e %12.6e',xqd(:,ip));
	 dxqd(:,ip)=ocf(tqd(ip),xt,ut,zt,upar);
    fprintf(ferr,'\n dx/dt %12.6e %12.6e %12.6e %12.6e %12.6e',dxqd(:,ip));
    for ig=1:ngc+1
     fprintf(ferr,'\n costates %12.6e %12.6e %12.6e %12.6e %12.6e',...
      psiqd((ig-1)*ns+1:ig*ns),ip);
    end
   end
   fprintf(ferr,'\n*********************************************');
   for ip=iqq-1:-1:1
    fprintf(ferr,'\n quadrature point no. %4i time =%10.5f',ip,tqd(ip));
    fprintf(ferr,'\n controls  %12.6e %12.6e %12.6e %12.6e %12.6e',uqd(:,ip));
    fprintf(ferr,'\n states  %12.6e %12.6e %12.6e %12.6e %12.6e',xqd(:,ip));
    fprintf(ferr,'\n dx/dt %12.6e %12.6e %12.6e %12.6e %12.6e',dxqd(:,ip));
    for ig=1:ngc+1
     fprintf(ferr,'\n costates %12.6e %12.6e %12.6e %12.6e %12.6e',...
      psiqd((ig-1)*ns+1:ig*ns),ip);
    end
   end
% stop execution
   fclose('all');
   error('Stopping Execution of MISER - all open files have been closed');
		 
  else
   psiqd(:,kpright:-1:kpleft)=yde';
   psit=psiqd(:,kpleft);
  end
 end
 lccostate=1;
end    

%_____________________________________________________________________________
% Compute quadratures of del H(i) del u, i=0,1,2,...,ngc.
if lgradf0
 gradf0=zeros(1,nctot);
 gf0error=zeros(1,nctot);
end
if lgradg
 gradg=zeros(ngc,nctot);
 ggerror=zeros(ngc,nctot);
end

if lgradf0 | lgradg
	
% Compute terms del x0 del z and del phi del z.
 if nz > 0
  if ~numder(3,5)
   dx0dz=ocdx0dz(zt,upar);
  else
   dx0dz=numdx0(zt,upar,reltest);
  end
  if lgradf0
   if ~numder(3,3)
    dpdz=ocdpdz(0,tfinal,xqd(:,nqp),zt,upar,rhoab,kabs);
   else
    dpdz=numdp(3,0,tfinal,xqd(:,nqp),zt,upar,rhoab,kabs,reltest);
   end	
   gradf0(nct+1:nctot)=dpdz+psiqd(1:ns,1)'*dx0dz;
  end
  if lgradg
   for ig=1:ngc
    if ~numder(3,3)
     dpdz=ocdpdz(ig,tau(ig),xqd(:,kpqtau(ig)),zt,upar,rhoab,kabs);
    else
     dpdz=numdp(3,ig,tau(ig),xqd(:,kpqtau(ig)),zt,upar,rhoab,kabs,reltest);
    end	
    gradg(ig,nct+1:nctot)=dpdz + psiqd(ig*ns+1:ns*(ig+1),1)'*dx0dz;
   end
  end
 end

%_____________________________________________________________________________
% Quadratures, one sub-interval at a time.
 if nc > 0
  ind=max(kpqk(2:nknot)-kpqk(1:nknot-1));
  if lgradf0, dh0duj=zeros(nct+1,ind+1); end
  if lgradg, dhiduj=zeros(ngc,nct+1,ind+1); end
 end	
 for ik=1:nknot-1
  kint=ik;
  dt=dtk(ik);
  numints=kpqk(ik+1)-kpqk(ik);
% Compute Romberg weights
  w=[3.5 16 6 16];
  werror=[-0.25 1 -1.5 1];
  for i=2:floor(numints/4);
   w=[w 7 16 6 16];
   werror=[werror -0.5 1 -1.5 1];
  end
  w=4*dt*[w 3.5]/45;
  werror=4*dt*[werror -0.25]/45;
  for iqq=1:numints+1
   iq=kpqk(ik)+iqq-1;
   tt=tqd(iq);
   xt=xqd(:,iq);
	
   if nc > 0
    ut = occalu(tt,kint,kpar,koncts,tk,kptk,kpcp);
    if ~numder(2,1)
     dfdu = ocdfdu(tt,xt,ut,zt,upar);
    else
     dfdu = numdf(2,tt,xt,ut,zt,upar,reltest);
    end
    if lgradf0
     if ~numder(2,2)
      df0du = ocdg0du(tt,xt,ut,zt,upar,rhoab,kabs);
     else
      df0du = numdg0(2,tt,xt,ut,zt,upar,rhoab,kabs,reltest);
     end
    end			
    if lgradg
     if ~numder(2,1)
      dgdu = ocdgdu(tt,xt,ut,zt,upar,rhoab,kabs);
     else
      dgdu = numdg(2,tt,xt,ut,zt,upar,rhoab,kabs,reltest);
     end
% For epsilon - tau algorithm.
     dgedg = ocg(tt,xt,ut,zt,upar,rhoab,kabs);
     for ig=1:ngc
      if lget(ig)
       dgedg(ig)=ocdetsm(dgedg(ig),epsjt);
      else
       dgedg(ig)=1;
      end
     end
    end			
    for ic=1:nc
% Compute integrands for derivatives of objective and constraints	wrt uj
     if lgradf0
      r0 = df0du(ic)+psiqd(1:ns,iq)'*dfdu(:,ic);
     end
     if lgradg
      for ig=1:ngc
       if kpqk(ik+1) > kpqtau(ig)
        r(ig)=0;
       else
        r(ig)=dgdu(ig,ic)*dgedg(ig)+psiqd(ig*ns+1:ns*(ig+1),iq)'*dfdu(:,ic);
       end
      end
     end
% Piecewise constant basis functions.
     inct=kpar(ic,ik);
     if koncts(ic) == 0
      if lgradf0, dh0duj(inct,iqq) = r0; end
      if lgradg, dhiduj(1:ngc,inct,iqq) = r; end							
% Piecewise linear continuous control.
     elseif koncts(ic) == 1
      tkleft=tk(kptk(ic)+inct-kpcp(ic));
      tkright=tk(kptk(ic)+inct-kpcp(ic)+1);
      brt=(tt-tkright)/(tkleft-tkright);
      blt=(tt-tkleft)/(tkright-tkleft);
      if lgradf0
       dh0duj(inct,iqq) = r0*brt;
       dh0duj(inct+1,iqq) = r0*blt;
      end						
      if lgradg
       for ig=1:ngc
        dhiduj(ig,inct,iqq) = r(ig)*brt;
        dhiduj(ig,inct+1,iqq) = r(ig)*blt;
       end
      end
     end
    end	% of ic loop
   end	% of nc > 0

%_____________________________________________________________________________
% Compute gradients with respect to system parameters.
   if nz > 0
    if ~numder(3,1)
     dfdz = ocdfdz(tt,xt,ut,zt,upar);
    else
     dfdz = numdf(3,tt,xt,ut,zt,upar,reltest);
    end
    if lgradf0
     if ~numder(3,2)
      df0dz = ocdg0dz(tt,xt,ut,zt,upar,rhoab,kabs);
     else
      df0dz = numdg0(3,tt,xt,ut,zt,upar,rhoab,kabs,reltest);
     end
    end			
    if lgradg
     if ~numder(3,4)
      dgdz = ocdgdz(tt,xt,ut,zt,upar,rhoab,kabs);
     else
      dgdz = numdg(3,tt,xt,ut,zt,upar,rhoab,kabs,reltest);
     end
% For epsilon - tau algorithm.
     dgedg = ocg(tt,xt,ut,zt,upar,rhoab,kabs);
     for ig=1:ngc
      if lget(ig)
       dgedg(ig)=ocdetsm(dgedg(ig),epsjt);
      else
       dgedg(ig)=1;
      end
     end
    end
% Compute integrands for derivatives of objective and constraints	wrt z
    if lgradf0, dh0dz(1:nz,iqq) = df0dz'+dfdz'*psiqd(1:ns,iq); end
    if lgradg
     for ig=1:ngc
      if kpqk(ik+1)>kpqtau(ig)
       dhidz(ig,1:nz,iqq)=zeros(1,nz,1);
      else
       dhidz(ig,1:nz,iqq)=dgdz(ig,1:nz)'*dgedg(ig) + ...
        dfdz'*psiqd(ig*ns+1:ns*(ig+1),iq);
      end
     end
    end
   end % of nz > 0
			
  end	% of iqq loop   
		
%_____________________________________________________________________________
% Compute Romberg quadrature and approx error for this interval
% control parameters first.
  for ic=1:nc
   inct=kpar(ic,ik);
   if lgradf0
    gradf0(inct)=gradf0(inct)+w*dh0duj(inct,1:numints+1)';
    gf0error(inct)=gf0error(inct)+werror*dh0duj(inct,1:numints+1)';
   end
   if lgradg
    for ig=1:ngc
     temp1=reshape(dhiduj(ig,inct,1:numints+1),1,numints+1);
     gradg(ig,inct)=gradg(ig,inct)+w*temp1';
     ggerror(ig,inct)=ggerror(ig,inct)+werror*temp1';
    end					
   end
   if koncts(ic)==1
    inct=inct+1;
    if lgradf0
     gradf0(inct)=gradf0(inct)+w*dh0duj(inct,1:numints+1)';
     gf0error(inct)=gf0error(inct)+werror*dh0duj(inct,1:numints+1)';
    end
    if lgradg
     for ig=1:ngc
      temp1=reshape(dhiduj(ig,inct,1:numints+1),1,numints+1);
      gradg(ig,inct)=gradg(ig,inct)+w*temp1';
      ggerror(ig,inct)=ggerror(ig,inct)+werror*temp1';
     end					
    end
   end
  end
%_____________________________________________________________________________
% Quadratures for system parameters.
  if nz > 0
   if lgradf0
    gradf0(nct+1:nctot)=gradf0(nct+1:nctot)+w*dh0dz(1:nz,1:numints+1)';
    gf0error(nct+1:nctot)=gf0error(nct+1:nctot)+werror*dh0dz(1:nz,1:numints+1)';
   end
   if lgradg
    for ig=1:ngc
     temp3=reshape(dhidz(ig,1:nz,1:numints+1),nz,numints+1);
     gradg(ig,nct+1:nctot)=gradg(ig,nct+1:nctot)+w*temp3';
     ggerror(ig,nct+1:nctot)=ggerror(ig,nct+1:nctot)+werror*temp3';
    end			
   end
  end

 end	% of for loop in ik

 lcgradf0=max([lcgradf0,lgradf0]);
 lcgradg=max([lcgradg,lgradg]);
		
end	% of (if lgradf0 | lgradg)

%_____________________________________________________________________________
% Set values into formal parameters
if lf0, f0=fj; end
if lg, g=fg; end
if lgradf0, gf0=gradf0; end
if lgradg, gg=gradg; end

Ntimes=Ntimes+1;

% End of ocfgeval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&

function u = occalu(t,kint,kpar,koncts,tk,kptk,kpcp)
% Calculate u at time t.
global Cb

nc=length(koncts);
for ic=1:nc
	icp=kpar(ic,kint);
	if koncts(ic)==0
		u(ic)=Cb(icp);
	else
		leftindex=kptk(ic)+icp-kpcp(ic);
		tkleft=tk(leftindex);
		tkright=tk(leftindex+1);
		u(ic)=Cb(icp)+(Cb(icp+1)-Cb(icp))*(t-tkleft)/(tkright-tkleft);
 	end
end
u=u';

% End of occalu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&

function [psi,lcint] = ocpsiset(x,z,kpqtf,ngc,kpright,kpqtau,tau,tfinal,...
psit,rhoab,kabs,upar,numder,reltest)	
% Computes initial values of costate equations and sets interval indicator to
% left or right of tau(i) values.

ns=length(x);
% Objective function costates
% kpqtf is pointer to last knot
if kpright == kpqtf
 if ~numder(1,3)
  dpdx = ocdpdx(0,tfinal,x,z,upar,rhoab,kabs);
 else
  dpdx = numdp(1,0,tfinal,x,z,upar,rhoab,kabs,reltest);
 end
 psi=dpdx';
 lcint(1)=1;
elseif kpright < kpqtf
 lcint(1)=1;
 psi=psit(1:ns,1);
else
 psi=zeros(ns,1); 
 lcint(1)=0;
end

% Constraint function costates
if ngc > 0
 for ig=1:ngc
  if kpright==kpqtau(ig)
   if ~numder(1,3)
    dpdx = ocdpdx(ig,tau(ig),x,z,upar,rhoab,kabs);
   else
    dpdx = numdp(1,ig,tau(ig),x,z,upar,rhoab,kabs,reltest);
   end
   psi = [psi; dpdx'];
   lcint(ig+1)=1;
  elseif kpright < kpqtau(ig)
   lcint(ig+1)=1;
   psi = [psi; psit(ns*ig+1:ns*(ig+1),1)];
  else
   psi=[psi; zeros(ns,1)];
   lcint(ig+1)=0;
  end
 end
end

% End of ocpsiset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&

function geps = ocetsm(g,epsjt)
% The epsilon-tau smoothing of min(0,g) as in Teo and Jennings 1989,
% for >= constraints.

if g >= epsjt
 geps=0;
elseif g <= -epsjt
 geps=g;
else
 geps=-0.25*(g-epsjt)^2/epsjt;
end

% End of ocetsm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&

function [worst,kflag,oworst,okflag,indf0,indf,indg,indp,indx0] = ...
 ocdcheck(t,x,u,z,worst,kflag,oworst,okflag,indf0,indf,indg,indp,...
 indx0,ngc,tstart,tfinal,tau,rhoab,kabs,reltest,ferr,kout,upar,numder)
% This routine computes second order numerical approximations to all
%  partial derivatives for which analytic derivatives are supplied and 
%  compares to the computed analytic ones at the current time.
%  The derivatives of phi and x0 are tested if t is appropriate.
%  Array kflag is returned, 0 - OK, 1 - fuzzy, 2 - wrong user derivatives.
%  okflag is the worst of kflag entries.
%  Array worst is returned with worst partial derivative.
%  oworst is the worst of the worst array.

ngcp=ngc+1;
ha=eps^(1/3);
ns=length(x);
nc=length(u);
nz=length(z);

% Derivatives wrt states
if any(~numder(1,1:4))
 if ~numder(1,2), df0 = ocdg0dx(t,x,u,z,upar,rhoab,kabs); end
 if ~numder(1,1), df = ocdfdx(t,x,u,z,upar); end
 if ngc > 0 & ~numder(1,4)
  dg = ocdgdx(t,x,u,z,upar,rhoab,kabs);
  gmax=zeros(ngc,1);
 end
 if ~numder(1,3)
  if t == tfinal, dp(1,:)=ocdpdx(0,t,x,z,upar,rhoab,kabs); end
  for ig=1:ngc
   if t == tau(ig), dp(ig+1,:)=ocdpdx(ig,t,x,z,upar,rhoab,kabs); end
  end
 end
 f0max=0;
 fmax=zeros(ns,1);
 pmax=zeros(1,ngc+1);

 for is=1:ns
  hr=ha*max([reltest abs(x(is))]);
  honhr=1/(2*hr);
  vsave=x(is);
% Forward function values.
  x(is)=vsave+hr;
  if ~numder(1,2), f0ph = ocg0(t,x,u,z,upar,rhoab,kabs); end
  if ~numder(1,1), fph = ocf(t,x,u,z,upar); end
  if ~numder(1,3) & t == tfinal
   pph(1) = ocphi(0,t,x,z,upar,rhoab,kabs);
  end
  if ngc > 0
   if ~numder(1,4), gph = ocg(t,x,u,z,upar,rhoab,kabs); end
   if ~numder(1,3)
    for ig=1:ngc
     if t == tau(ig), pph(ig+1) = ocphi(ig,t,x,z,upar,rhoab,kabs); end
    end
   end
  end
% Backward function values and derivatives
  x(is)=vsave-hr;
% Derivatives of f0 wrt state x.
  if ~numder(1,2)
   f0mh = ocg0(t,x,u,z,upar,rhoab,kabs);
   f0max=max([f0max abs(f0ph) abs(f0mh)]);
   df0(is)=df0(is)-(f0ph-f0mh)*honhr;
  end
% Derivatives of f wrt state x.
  if ~numder(1,1)
   fmh = ocf(t,x,u,z,upar);
   fmax=(max([fmax';abs(fph');abs(fmh')]))';
   df(:,is)=df(:,is)-(fph-fmh)*honhr;
  end
% Derivatives of g wrt state x.
  if ngc > 0 & ~numder(1,4)
   gmh = ocg(t,x,u,z,upar,rhoab,kabs);
   gmax=(max([gmax';abs(gph');abs(gmh')]))';
   dg(:,is)=dg(:,is)-(gph-gmh)*honhr;
  end
% Derivative of phi wrt state x.
  if ~numder(1,3)
   if t == tfinal
    pmh(1) = ocphi(0,t,x,z,upar,rhoab,kabs);
    pmax(1)=max([pmax(1) abs(pph(1)),abs(pmh(1))]);
    dp(1,is)=dp(1,is)-(pph(1)-pmh(1))*honhr;
   end
   for ig=1:ngc
    if t == tau(ig)
     pmh(ig+1) = ocphi(ig,t,x,z,upar,rhoab,kabs);
     pmax(ig+1)=max([pmax(ig+1) abs(pph(ig+1)) abs(pmh(ig+1))]);
     dp(ig+1,is)=dp(ig+1,is)-(pph(ig+1)-pmh(ig+1))*honhr;
    end
   end
  end
  x(is)=vsave;
 end % of is loop

%_____________________________________________________________________________
% Check derivative of f0 wrt state x.
 if ~numder(1,2)
  for is=1:ns
   [worst(1,1) kflag(1,1) indf0(is)] = ocdprint(t,x,u,z,df0(is),...
    worst(1,1),kflag(1,1),'objective ',0,'state  ',is,f0max,...
    indf0(is),reltest,ferr,kout);
  end
 end
% Check derivatives of f wrt state x.
 if ~numder(1,1)
  for ifs=1:ns
   for is=1:ns
    [worst(2,1) kflag(2,1) indf(ifs,is)] = ocdprint(t,x,u,z,...
     df(ifs,is),worst(2,1),kflag(2,1),'state     ',ifs,'state  ',...
     is,fmax(ifs),indf(ifs,is),reltest,ferr,kout);
   end
  end
 end
% Check derivatives of g wrt state x.
 if ~numder(1,4)
  for ig=1:ngc
   for is=1:ns
    [worst(3,1) kflag(3,1) indg(ig,is)] = ocdprint (t,x,u,z,dg(ig,is),...
     worst(3,1),kflag(3,1),'constraint',ig,'state  ',is,gmax(ig),...
     indg(ig,is),reltest,ferr,kout);
   end
  end
 end
% Check derivatives of phi wrt state x.
 if ~numder(1,3)
  if t == tfinal
   for is=1:ns
    [worst(4,1) kflag(4,1) indp(1,is)] = ocdprint(t,x,u,z,dp(1,is),...
     worst(4,1),kflag(4,1),'phi      ',0,'state  ',is,pmax(1),...
     indp(1,is),reltest,ferr,kout);
   end
  end
  for ig=1:ngc
   if t == tau(ig)
    for is=1:ns
     [worst(4,1) kflag(4,1) indp(ig+1,is)] = ocdprint(t,x,u,z,...
      dp(ig+1,is),worst(4,1),kflag(4,1),'phi      ',ig,'state  ',...
      is,pmax(ig+1),indp(ig+1,is),reltest,ferr,kout);
    end
   end
  end
 end

end % of derivatives wrt states
%_____________________________________________________________________________
% Derivatives wrt controls.
if nc > 0 & (numder(2,1)==0 | numder(2,2)==0 | numder(2,4)==0)
 if ~numder(2,2), df0 = ocdg0du(t,x,u,z,upar,rhoab,kabs); end
 if ~numder(2,1), df = ocdfdu(t,x,u,z,upar); end
 if ngc > 0 & ~numder(2,4)
  dg = ocdgdu(t,x,u,z,upar,rhoab,kabs);
  gmax=zeros(ngc,1);
 end
 f0max=0;
 fmax=zeros(ns,1);
 for ic=1:nc
  hr=ha*max([reltest abs(u(ic))]);
  honhr=1/(hr+hr);
  vsave=u(ic);
% Forward function values.
  u(ic)=vsave+hr;
  if ~numder(2,2), f0ph = ocg0(t,x,u,z,upar,rhoab,kabs); end
  if ~numder(2,1), fph = ocf(t,x,u,z,upar); end
  if ngc > 0 & ~numder(2,4) , gph = ocg(t,x,u,z,upar,rhoab,kabs); end
% Backward function values and derivatives
  u(ic)=vsave-hr;
% Derivatives of f0 wrt control u.
  if ~numder(2,2)
   f0mh = ocg0(t,x,u,z,upar,rhoab,kabs);
   f0max=max([f0max abs(f0ph) abs(f0mh)]);
   df0(ic)=df0(ic)-(f0ph-f0mh)*honhr;
  end
% Derivatives of f wrt control u.
  if ~numder(2,1)
   fmh = ocf(t,x,u,z,upar);
   fmax=(max([fmax';abs(fph');abs(fmh')]))';
   df(:,ic)=df(:,ic)-(fph-fmh)*honhr;
  end
% Derivatives of g wrt control u.
  if ngc > 0 & ~numder(2,4)
   gmh = ocg(t,x,u,z,upar,rhoab,kabs);
   gmax=(max([gmax';abs(gph');abs(gmh')]))';
   dg(:,ic)=dg(:,ic)-(gph-gmh)*honhr;
  end
  u(ic)=vsave;
 end
	
%_____________________________________________________________________________
% Check derivative of f0 wrt control u.
 if ~numder(2,2)
  for ic=1:nc
   [worst(1,2) kflag(1,2) indf0(ns+ic)] = ocdprint(t,x,u,z,df0(ic),...
    worst(1,2),kflag(1,2),'objective ',0,'control',ic,f0max,...
    indf0(ns+ic),reltest,ferr,kout);
  end
 end
% Check derivatives of f wrt control u.
 if ~numder(2,1)
  for ifc=1:ns
   for ic=1:nc
    [worst(2,2) kflag(2,2) indf(ifc,ns+ic)] = ocdprint(t,x,u,z,df(ifc,ic),...
     worst(2,2),kflag(2,2),'state     ',ifc,'control',ic,fmax(ifc),...
     indf(ifc,ns+ic),reltest,ferr,kout);
   end
  end
 end
% Check derivatives of g wrt control u.
 if ~numder(2,4)
  for ig=1:ngc
   for ic=1:nc
    [worst(3,2) kflag(3,2) indg(ig,ns+ic)] = ocdprint (t,x,u,z,dg(ig,ic),...
     worst(3,2),kflag(3,2),'constraint',ig,'control',ic,gmax(ig),...
     indg(ig,ns+ic),reltest,ferr,kout);
   end
  end
 end
	
end

%_____________________________________________________________________________
% Derivatives wrt system parameters.
if nz > 0 & any(~numder(3,:))
 if ~numder(3,2), df0 = ocdg0dz(t,x,u,z,upar,rhoab,kabs); end
 if ~numder(3,1), df = ocdfdz(t,x,u,z,upar); end
 if ngc > 0 & ~numder(3,4)
  dg = ocdgdz(t,x,u,z,upar,rhoab,kabs);
  gmax=zeros(ngc,1);
 end
 if t == tstart & ~numder(3,5), dx0=ocdx0dz(z,upar); end
 if ~numder(3,3)
  if t == tfinal, dp(1,1:nz)=ocdpdz(0,t,x,z,upar,rhoab,kabs); end
  for ig=1:ngc
   if t == tau(ig), dp(ig+1,1:nz)=ocdpdz(ig,t,x,z,upar,rhoab,kabs); end
  end
 end
 f0max=0;
 fmax=zeros(ns,1);
 pmax=zeros(1,ngc+1);
 x0max=zeros(ns,1);
 
 for iz=1:nz
  hr=ha*max([reltest abs(z(iz))]);
  honhr=1/(hr+hr);
  vsave=z(iz);
% Forward function values.
  z(iz)=vsave+hr;
  if ~numder(3,2), f0ph = ocg0(t,x,u,z,upar,rhoab,kabs); end
  if ~numder(3,1), fph = ocf(t,x,u,z,upar); end
  if t == tstart & ~numder(3,5), x0ph=ocxzero(z,upar); end
  if t == tfinal & ~numder(3,3), pph(1) = ocphi(0,t,x,z,upar,rhoab,kabs); end
  if ngc > 0
   if ~numder(3,4), gph = ocg(t,x,u,z,upar,rhoab,kabs); end
   if ~numder(3,3)
    for ig=1:ngc
     if t == tau(ig), pph(ig+1) = ocphi(ig,t,x,z,upar,rhoab,kabs); end
    end
   end
  end
% Backward function values and derivatives
  z(iz)=vsave-hr;
% Derivatives of f0 wrt system parameters z.
  if ~numder(3,2)
   f0mh = ocg0(t,x,u,z,upar,rhoab,kabs);
   f0max=max([f0max abs(f0ph) abs(f0mh)]);
   df0(iz)=df0(iz)-(f0ph-f0mh)*honhr;
  end
% Derivatives of f wrt system parameters.
  if ~numder(3,1)
   fmh = ocf(t,x,u,z,upar);
   fmax=(max([fmax';abs(fph');abs(fmh')]))';
   df(:,iz)=df(:,iz)-(fph-fmh)*honhr;
  end
% Derivatives of g wrt system parameters z.
  if ngc > 0 & ~numder(3,4)
   gmh = ocg(t,x,u,z,upar,rhoab,kabs);
   gmax=(max([gmax';abs(gph');abs(gmh')]))';
   dg(:,iz)=dg(:,iz)-(gph-gmh)*honhr;
  end
% Derivative of phi wrt system parameters z.
  if t == tfinal & ~numder(3,3)
   pmh(1) = ocphi(0,t,x,z,upar,rhoab,kabs);
   pmax(1)=max([pmax(1) abs(pph(1)) abs(pmh(1))]);
   dp(1,iz)=dp(1,iz)-(pph(1)-pmh(1))*honhr;
  end
  if ~numder(3,3)
   for ig=1:ngc
    if t == tau(ig)
     pmh(ig+1) = ocphi(ig,t,x,z,upar,rhoab,kabs);
     pmax(ig+1)=max([pmax(ig+1) abs(pph(ig+1)) abs(pmh(ig+1))]);
     dp(ig+1,iz)=dp(ig+1,iz)-(pph(ig+1)-pmh(ig+1))*honhr;
    end
   end
  end
% Derivative of x0 wrt system parameters z.
  if t == tstart & ~numder(3,5)
   x0mh= ocxzero(z,upar);
   x0max=(max([x0max';abs(x0ph');abs(x0mh')]))';
   dx0(:,iz)=dx0(:,iz)-(x0ph-x0mh)*honhr;
  end
  z(iz)=vsave;
 end % of iz loop
	
%_____________________________________________________________________________
% Check derivative of f0 wrt system parameters z.
 if ~numder(3,2)
  for iz=1:nz
   [worst(1,3) kflag(1,3) indf0(ns+nc+iz)] = ocdprint(t,x,u,z,...
    df0(iz),worst(1,3),kflag(1,3),'objective ',0,'system ',iz,...
    f0max,indf0(ns+nc+iz),reltest,ferr,kout);
  end
 end
% Check derivatives of f wrt system parameters z.
 if ~numder(3,1)
  for ifs=1:ns
   for iz=1:nz
    [worst(2,3) kflag(2,3) indf(ifs,ns+nc+iz)] = ocdprint(t,x,u,z,...
     df(ifs,iz),worst(2,3),kflag(2,3),'state     ',ifs,'system ',...
     iz,fmax(ifs),indf(ifs,ns+nc+iz),reltest,ferr,kout);
   end
  end
 end
% Check derivatives of g wrt system parameters z.
 if ~numder(3,4)
  for ig=1:ngc
   for iz=1:nz
    [worst(3,3) kflag(3,3) indg(ig,ns+nc+iz)] = ocdprint(t,x,u,z,...
     dg(ig,iz),worst(3,3),kflag(3,3),'constraint',ig,'system ',...
     iz,gmax(ig),indg(ig,ns+nc+iz),reltest,ferr,kout);
   end
  end
 end
% Check derivatives of phi wrt system parameters z.
 if ~numder(3,3)
  if t == tfinal
   for iz=1:nz
    [worst(4,3) kflag(4,3) indp(1,ns+iz)] = ocdprint(t,x,u,z,...
     dp(1,iz),worst(4,3),kflag(4,3),'phi      ',0,'system ',...
     iz,pmax(1),indp(1,ns+iz),reltest,ferr,kout);
   end
  end
  for ig=1:ngc
   if t == tau(ig)
    for iz=1:nz
     [worst(4,3) kflag(4,3) indp(ig+1,ns+iz)] = ocdprint(t,x,u,z,...
      dp(ig+1,iz),worst(4,3),kflag(4,3),'phi      ',ig,'system ',...
      iz,pmax(ig+1),indp(ig+1,ns+iz),reltest,ferr,kout);
     end
   end
  end
 end
% Check derivatives of x0 wrt system parameters z.
 if t == tstart & ~numder(3,5)
  for is=1:ns
   for iz=1:nz
    [worst(5,3) kflag(5,3) indx0(is,iz)] = ocdprint(t,x,u,z,...
     dx0(is,iz),worst(5,3),kflag(5,3),'  x0      ',is,'system ',...
     iz,x0max(is),indx0(is,iz),reltest,ferr,kout);
   end
  end
 end
	
end

%_____________________________________________________________________________
oworst=max([oworst max(worst)]);
okflag=max([okflag max(kflag)]);

% End of ocdcheck
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&

function [wor,kfl,inderr] = ocdprint(t,x,u,z,diff,wor,kfl,namef,ii,namev,...
 iv,fmax,inderr,reltest,ferr,kout)
%  Prints the error messages for ocdcheck.

ns=length(x);
nc=length(u);
nz=length(z);
ha=eps^(1/3);
srepsm=sqrt(eps);
rmax=max([reltest fmax]);
adiff=abs(diff);
wor=max([wor adiff]);
if adiff > ha*rmax
 fprintf(ferr,'\n OCDCHECK --- At t=%12.3e',t);
 fprintf(ferr,'\n partial derivative and centred difference DISAGREE  ');
 fprintf(ferr,'\n in %s %3i function, wrt %s variable number %3i',...
  namef,ii,namev,iv);
 fprintf(ferr,'\n partial - numerical = %15.7e',diff);
 kfl=2;
 inderr=inderr+1000;
 if kout(1) >= 2
  fprintf(ferr,'\n States:%12.4e%12.4e%12.4e%12.4e%12.4e',x);
  if nc > 0, fprintf(ferr,'\n Controls:%12.4e%12.4e%12.4e%12.4e%12.4e',u); end
  if nz > 0, fprintf(ferr,'\n System:%12.4e%12.4e%12.4e%12.4e%12.4e',z); end
 end
	
elseif adiff > srepsm*rmax
 if kout(1) >= 1
  fprintf(ferr,'\n OCDCHECK --- At t=%12.3e',t);
  fprintf(ferr,'\n partial derivative and centred difference is fuzzy  ');
  fprintf(ferr,'\n in %s %3i function, wrt %s variable number %3i',...
   namef,ii,namev,iv);
  fprintf(ferr,'\n partial - numerical = %15.7e',diff);
  kfl=max([kfl 1]);
  inderr=inderr+1;
 end
else
 if kout(1) >= 3
  fprintf(ferr,'\n OCDCHECK --- At t=%12.3e',t);
  fprintf(ferr,'\n partial derivative and centred difference  is OK    ');
  fprintf(ferr,'\n in %s %3i function, wrt %s variable number %3i',...
   namef,ii,namev,iv);
  fprintf(ferr,'\n partial - numerical = %15.7e',diff);
 end
end

% End of ocdprint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&

function ociprint(ns,nc,nz,nqp,indf0,indf,indg,indp,indx0,kflag,ferr)
% Prints pattern of errors in user derivatives from indicator arrays.

fprintf('\n Entries indicate number of points where errors are');
fprintf(' correct, fuzzy and wrong.\n The groups of derivatives');
fprintf(' shown are only those where the OK flag is non-zero.');
fprintf(ferr,'\n Entries indicate number of points where errors are');
fprintf(ferr,' correct, fuzzy and wrong.\n The groups of derivatives');
fprintf(ferr,' shown are only those where the OK flag is non-zero.');
nall=ns+nc+nz;
		
if nall <= 12
 if sum(kflag(1,1:3)) > 0
  ocjprint('objective ','      all ',indf0,ns,nc,nz,nqp,ferr);
 end
 if sum(kflag(2,1:3)) > 0
  ocjprint('    state ','      all ',indf,ns,nc,nz,nqp,ferr);
 end
 if sum(kflag(3,1:3)) > 0
  ocjprint('constraint','      all ',indg,ns,nc,nz,nqp,ferr);
 end
 if kflag(4,1) > 0
  ocjprint('      phi ','    state ',indp(:,1:ns),ns,nc,nz,nqp,ferr);
 end
 if kflag(4,3) > 0
  ocjprint('      phi ','   system ',indp(:,ns+1:ns+nz),ns,nc,nz,nqp,ferr);
 end
 if kflag(5,3) > 0
  ocjprint('      x0  ','   system ',indx0,ns,nc,nz,nqp,ferr);
 end
else
 if kflag(1,1)>0
  ocjprint('objective ','    state ',indf0(1:ns),ns,nc,nz,nqp,ferr);
 end
 if kflag(1,2)>0
  ocjprint('objective ','  control ',indf0(ns+1:ns+nc),ns,nc,nz,nqp,ferr);
 end
 if kflag(1,3)>0
  ocjprint('objective ','   system ',indf0(ns+nc+1:nall),ns,nc,nz,nqp,ferr);
 end
 if kflag(2,1)>0
  ocjprint('    state ','    state ',indf(:,1:ns),ns,nc,nz,nqp,ferr);
 end
 if kflag(2,2)>0
  ocjprint('    state ','  control ',indf(:,ns+1:ns+nc),ns,nc,nz,nqp,ferr);
 end
 if kflag(2,3)>0
  ocjprint('    state ','   system ',indf(:,ns+nc+1:nall),ns,nc,nz,nqp,ferr);
 end
 if kflag(3,1)>0
  ocjprint('constraint','    state ',indg(:,1:ns),ns,nc,nz,nqp,ferr);
 end
 if kflag(3,2)>0
  ocjprint('constraint','  control ',indg(:,ns+1:ns+nc),ns,nc,nz,nqp,ferr);
 end
 if kflag(3,3)>0
  ocjprint('constraint','   system ',indg(:,ns+nc+1:nall),ns,nc,nz,nqp,ferr);
 end	
 if kflag(4,1)>0
  ocjprint('      phi ','    state ',indp(:,1:ns),ns,nc,nz,nqp,ferr);
 end
 if kflag(4,3)>0
  ocjprint('      phi ','   system ',indp(:,ns+1:ns+nz),ns,nc,nz,nqp,ferr);
 end
 if kflag(5,3)>0
  ocjprint('      x0  ','   system ',indx0,ns,nc,nz,nqp,ferr);
 end
end
	
% End of ociprint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&

function ocjprint(side,top,ind,ns,nc,nz,nqp,ferr)

[numv numh]=size(ind);
if top == '      all '
 fprintf('\n\n                              all = state, control, system');
 fprintf(ferr,'\n\n                              all = state, control, system');
 fprintf('\n                   %5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i',...
  [1:ns 1:nc 1:nz]);
 fprintf(ferr,'\n                   %5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i',...
  [1:ns 1:nc 1:nz]); 
else
 fprintf('\n\n                    %s',top);
 fprintf(ferr,'\n\n                    %s',top);
 fprintf('\n                   %5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i',1:numh);
 fprintf(ferr,'\n                   %5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i',...
  1:numh);
end
tot=nqp;
if (side == '      phi ') | (side == '      x0  '), tot=1; end
for i=1:numv
 in=i;
 if side == '      phi ', in=i-1; end
 row3=floor(ind(i,:)/1000);
 row2=ind(i,:)-1000*row3;
 row1=tot-row2-row3;
 fprintf('\n\n                OK  %5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i',row1);
 fprintf(ferr,'\n                OK  %5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i',...
  row1);
 fprintf('\n%s%4i FUZZY%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i',side,in,row2);
 fprintf(ferr,'\n%s%4i FUZZY%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i',...
  side,in,row2);			
 fprintf('\n               WRONG%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i',row3);
 fprintf(ferr,'\n               WRONG%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i',...
  row3);			
end

% End of ocjprint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&

function df = numdf(sw,t,x,u,z,upar,reltest);
% Finite difference code for derivatives of f. Uses central differences.
% Not recommended for very accurate work, and it is slow.
% sw = 1  	dfdx
% sw = 2  	dfdu
% sw = 3  	dfdz

ha=eps^(1/3);
if sw==1
 nx=length(x);
 for i=1:nx
  xi=x(i);
  hh=ha*max(reltest,abs(xi));
  x(i)=xi+hh;
  fph = ocf(t,x,u,z,upar);
  x(i)=xi-hh;
  fmh = ocf(t,x,u,z,upar);
  df(:,i)=(fph-fmh)/(2*hh);
  x(i)=xi;
 end
elseif sw==2
 nu=length(u);
 for i=1:nu
  ui=u(i);
  hh=ha*max(reltest,abs(ui));
  u(i)=ui+hh;
  fph = ocf(t,x,u,z,upar);
  u(i)=ui-hh;
  fmh = ocf(t,x,u,z,upar);
  df(:,i)=(fph-fmh)/(2*hh);
  u(i)=ui;
 end
else
 nz=length(z);
 for i=1:nz
  zi=z(i);
  hh=ha*max(reltest,abs(zi));
  z(i)=zi+hh;
  fph = ocf(t,x,u,z,upar);
  z(i)=zi-hh;
  fmh = ocf(t,x,u,z,upar);
  df(:,i)=(fph-fmh)/(2*hh);
  z(i)=zi;
 end
end

% End of numdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&

function dg0 = numdg0(sw,t,x,u,z,upar,rhoab,kabs,reltest)
% Finite difference code for derivatives of g0. Uses central differences.
% Not recommended for very accurate work, and it is slow.
% sw = 1  	dg0dx
% sw = 2  	dg0du
% sw = 3  	dg0dz

ha=eps^(1/3);
if sw==1
 nx=length(x);
 for i=1:nx
  xi=x(i);
  hh=ha*max(reltest,abs(xi));
  x(i)=xi+hh;
  g0ph = ocg0(t,x,u,z,upar,rhoab,kabs);
  x(i)=xi-hh;
  g0mh = ocg0(t,x,u,z,upar,rhoab,kabs);
  dg0(i)=(g0ph-g0mh)/(2*hh);
  x(i)=xi;
 end
elseif sw==2
 nu=length(u);
 for i=1:nu
  ui=u(i);
  hh=ha*max(reltest,abs(ui));
  u(i)=ui+hh;
  g0ph = ocg0(t,x,u,z,upar,rhoab,kabs);
  u(i)=ui-hh;
  g0mh = ocg0(t,x,u,z,upar,rhoab,kabs);
  dg0(i)=(g0ph-g0mh)/(2*hh);
  u(i)=ui;
 end
else
 nz=length(z);
 for i=1:nz
  zi=z(i);
  hh=ha*max(reltest,abs(zi));
  z(i)=zi+hh;
  g0ph = ocg0(t,x,u,z,upar,rhoab,kabs);
  z(i)=zi-hh;
  g0mh = ocg0(t,x,u,z,upar,rhoab,kabs);
  dg0(i)=(g0ph-g0mh)/(2*hh);
  z(i)=zi;
 end
end

% End of numdg0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&

function dphi = numdp(sw,ig,taut,x,z,upar,rhoab,kabs,reltest);
% Finite difference code for derivatives of phi. Uses central differences.
% Not recommended for very accurate work, and it is slow.
% sw = 1  	dphidx
% sw = 3  	dphidz

ha=eps^(1/3);
if sw==1
 nx=length(x);
 for i=1:nx
  xi=x(i);
  hh=ha*max(reltest,abs(xi));
  x(i)=xi+hh;
  phiph = ocphi(ig,taut,x,z,upar,rhoab,kabs);
  x(i)=xi-hh;
  phimh = ocphi(ig,taut,x,z,upar,rhoab,kabs);
  dphi(i)=(phiph-phimh)/(2*hh);
  x(i)=xi;
 end
else
 nz=length(z);
 for i=1:nz
  zi=z(i);
  hh=ha*max(reltest,abs(zi));
  z(i)=zi+hh;
  phiph = ocphi(ig,taut,x,z,upar,rhoab,kabs);
  z(i)=zi-hh;
  phimh = ocphi(ig,taut,x,z,upar,rhoab,kabs);
  dphi(i)=(phiph-phimh)/(2*hh);
  z(i)=zi;
 end
end

% End of numdp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&

function dg = numdg(sw,t,x,u,z,upar,rhoab,kabs,reltest);
% Finite difference code for derivatives of g. Uses central differences.
% Not recommended for very accurate work, and it is slow.
% sw = 1  	dgdx
% sw = 2  	dgdu
% sw = 3  	dgdz

ha=eps^(1/3);
if sw==1
 nx=length(x);
 for i=1:nx
  xi=x(i);
  hh=ha*max(reltest,abs(xi));
  x(i)=xi+hh;
  gph = ocg(t,x,u,z,upar,rhoab,kabs);
  x(i)=xi-hh;
  gmh = ocg(t,x,u,z,upar,rhoab,kabs);
  dg(:,i)=(gph-gmh)/(2*hh);
  x(i)=xi;
 end
elseif sw==2
 nu=length(u);
 for i=1:nu
  ui=u(i);
  hh=ha*max(reltest,abs(ui));
  u(i)=ui+hh;
  gph = ocg(t,x,u,z,upar,rhoab,kabs);
  u(i)=ui-hh;
  gmh = ocg(t,x,u,z,upar,rhoab,kabs);
  dg(:,i)=(gph-gmh)/(2*hh);
  u(i)=ui;
 end
else
 nz=length(z);
 for i=1:nz
  zi=z(i);
  hh=ha*max(reltest,abs(zi));
  z(i)=zi+hh;
  gph = ocg(t,x,u,z,upar,rhoab,kabs);
  z(i)=zi-hh;
  gmh = ocg(t,x,u,z,upar,rhoab,kabs);
  dg(:,i)=(gph-gmh)/(2*hh);
  z(i)=zi;
 end
end

% End of numdg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&

function dx0 = numdx0(z,upar,reltest);
% Finite difference code for derivative of x0 wrt z. Uses central differences.
% Not recommended for very accurate work, and it is slow.

ha=eps^(1/3);
nz=length(z);
for i=1:nz
 zi=z(i);
 hh=ha*max(reltest,abs(zi));
 z(i)=zi+hh;
 x0ph = ocxzero(z,upar);
 z(i)=zi-hh;
 x0mh = ocxzero(z,upar);
 dx0(:,i)=(x0ph-x0mh)/(2*hh);
 z(i)=zi;
end

% End of numdx0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&
