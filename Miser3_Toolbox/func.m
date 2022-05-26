function [f0,gradf0] = func(c,ocsame,tolx,reltest,kpqk,tknot,dtk,h,tk,...
 tqd,tau,kpqtau,kpar,koncts,kpcp,kncp,kptk,ns,nz,lget,numder,epsjt,taujt,...
 rhoab,tolpsi,tfinal,ferr,kout,consdata1,misc,kreg,wreg0,wreg1,wreg2,lreg);
% routine to compute objective and gradient for FMINCON

ngc=consdata1(3)+consdata1(6);
nct=length(c)-nz;
kabs=misc(5);
if nargout==1
 flag=0002;
 if ocsame, flag=0001; end
	
 f0 = ocfgeval(flag,c,tolx,reltest,kpqk,tknot,...
  dtk,tk,tqd,tau,kpqtau,kpar,koncts,kpcp,kptk,ns,nz,ngc,lget,numder,...
  epsjt,taujt,rhoab,tolpsi,tfinal,ferr,kout,misc);
		
 if lreg
  f0 = ocregzn(f0,c,kpcp,koncts,kncp,h,kptk,kreg,wreg0,wreg1,wreg2,rhoab,kabs);
 end
	
else
 flag=0202;
 if ocsame, flag=0101; end

 [f0,gradf0]= ocfgeval(flag,c,tolx,reltest,kpqk,tknot,...
  dtk,tk,tqd,tau,kpqtau,kpar,koncts,kpcp,kptk,ns,nz,ngc,lget,numder,...
  epsjt,taujt,rhoab,tolpsi,tfinal,ferr,kout,misc);

 if lreg 
  f0 = ocregzn(f0,c,kpcp,koncts,kncp,h,kptk,kreg,wreg0,wreg1,wreg2,rhoab,kabs);
  gradf0 = ocdregzn(gradf0,c,kpcp,koncts,kncp,h,kptk,kreg,wreg0,wreg1,...
   wreg2,rhoab,kabs);
 end
	
end

% End of func	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
		
function f = ocregzn(f,c,kpcp,koncts,kncp,h,kptk,kreg,wreg0,wreg1,wreg2,...
 rhoab,kabs)
% adds regularization terms to the objective

for ic=1:length(koncts)
 for k=1:3:7
  order=kreg(ic,k);
  if order~=-1
   norm=kreg(ic,k+1);
   kcpt=kpcp(ic)-1;
			
   if order==0
    if norm==1
     if koncts(ic)==0
      for j=1:kncp(ic)
       f=f+wreg0(ic,j)*ocabsm(c(kcpt+j),rhoab,kabs);
      end
     elseif koncts(ic)==1
      for j=2:kncp(ic)
       if c(kcpt+j-1)*c(kcpt+j)>=0
        f=f+wreg0(ic,j)*ocabsm(c(kcpt+j-1)+c(kcpt+j),rhoab,kabs);
       else
        f=f+wreg0(ic,j)*ocabsm((c(kcpt+j-1)^2+...
         c(kcpt+j)^2)/(c(kcpt+j)-c(kcpt+j-1)),rhoab,kabs);
       end
      end
     end
    elseif norm==2
     if koncts(ic)==0
      for j=1:kncp(ic)
       f=f+wreg0(ic,j)*c(kcpt+j)^2;
      end
     elseif koncts(ic)==1
      for j=2:kncp(ic)
       f=f+wreg0(ic,j)*(c(kcpt+j-1)^2+c(kcpt+j-1)*c(kcpt+j)+c(kcpt+j)^2);
      end
     end
    end
			
   elseif order==1
    if norm==1
     for j=2:kncp(ic)
      f=f+wreg1(ic,j)*ocabsm(c(kcpt+j)-c(kcpt+j-1),rhoab,kabs);
     end
    elseif norm==2
     for j=2:kncp(ic)
      f=f+wreg1(ic,j)*(c(kcpt+j)-c(kcpt+j-1))^2;
     end
    end
				
   elseif order==2
    hj=h(kptk(ic)+1);
    if norm==1
     for j=2:kncp(ic)-1
      hjp=h(kptk(ic)+j);
      vjp=(c(kcpt+j+1)-c(kcpt+j))/hjp-(c(kcpt+j)-c(kcpt+j-1))/hj;
      f=f+wreg2(ic,j)*ocabsm(vjp,rhoab,kabs);
      hj=hjp;
     end
    elseif norm==2
     for j=2:kncp(ic)-1
      hjp=h(kptk(ic)+j);
      vjp=(c(kcpt+j+1)-c(kcpt+j))/hjp-(c(kcpt+j)-c(kcpt+j-1))/hj;
      f=f+wreg2(ic,j)*vjp^2;
      hj=hjp
     end
    end
   end
  end
		
 end
end
		
% End of ocregzn

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function df = ocdregzn (df,c,kpcp,koncts,kncp,h,kptk,kreg,wreg0,...
 wreg1,wreg2,rhoab,kabs)
% adds gradient of regularization terms to the objective gradient

for ic=1:length(koncts)
 for k=1:3:7
  order=kreg(ic,k);
  if order~=-1
   norm=kreg(ic,k+1);
   kcpt=kpcp(ic)-1;
			
   if order==0
    if norm==1
     if koncts(ic)==0
      for j=1:kncp(ic)
       df(kcpt+j)=df(kcpt+j)+wreg0(ic,j)*ocdabsm(c(kcpt+j),rhoab,kabs);
      end
     elseif koncts(ic)==1
      sj=c(kcpt+1);
      for j=1:kncp(ic)
       if j~=kncp(ic)
        sjp=c(kcpt+j+1);
        wjp=wreg0(ic,j+1);
        if sj*sjp>=0
         df(kcpt+j)=df(kcpt+j)+wjp*ocdabsm(sj+sjp,rhoab,kabs);
        else
         df(kcpt+j)=df(kcpt+j)+wjp*ocdabsm((sj^2+sjp^2)/(sjp-sj),rhoab,kabs)...
          *(sjp^2-sj^2+2*sj*sjp)/(sjp-sj)^2;
        end
       end
       if j~=1
        if sj*sjm>=0
         df(kcpt+j)=df(kcpt+j)+wj*ocdabsm(sj+sjm,rhoab,kabs);
        else
         df(kcpt+j)=df(kcpt+j)+wj*ocdabsm((sj^2+sjm^2)/(sj-sjm),rhoab,kabs)...
          *(sj^2-sjm^2-2*sj*sjm)/(sj-sjm)^2;
        end
       end
       sjm=sj;
       sj=sjp;
       wj=wjp;
      end
     end
    elseif norm==2
     if koncts(ic)==0
      for j=1:kncp(ic)
       df(kcpt+j)=df(kcpt+j)+wreg0(ic,j)*2*(c(kcpt+j));
      end
     elseif koncts(ic)==1
      sj=c(kcpt+1);
      for j=1:kncp(ic)
       if j~=kncp(ic)
        sjp=c(kcpt+j+1);
        wj=wreg0(ic,j+1);
        df(kcpt+j)=df(kcpt+j)+wj*(2*sj+sjp);
       end
       if j~=1, df(kcpt+j)=df(kcpt+j)+wjm*(2*sj+sjm); end
       sjm=sj;
       sj=sjp;
       wjm=wj;
      end
     end
    end
				
   elseif order==1
    if norm==1
     sj=c(kcpt+1);
     wj=wreg1(ic,1);
     for j=1:kncp(ic)
      if j~=kncp(ic)
       sjp=c(kcpt+j+1);
       wjp=wreg1(ic,j+1);
       df(kcpt+j)=df(kcpt+j)-wjp*ocdabsm(sjp-sj,rhoab,kabs);
      end
      if j~=1, df(kcpt+j)=df(kcpt+j)+wj*ocdabsm(sj-sjm,rhoab,kabs); end
      sjm=sj;
      sj=sjp;
      wj=wjp;
     end
    elseif norm==2
     sj=c(kcpt+1);
     wj=wreg1(ic,1)*2;
     for j=1:kncp(ic)
      if j~=kncp(ic)
       sjp=c(kcpt+j+1);
       wjp=wreg1(ic,j+1);
       wjp=wreg1(ic,j+1)*2;
       df(kcpt+j)=df(kcpt+j)-wjp*(sjp-sj);
      end
      if j~=1, df(kcpt+j)=df(kcpt+j)+wj*(sj-sjm); end
      sjm=sj;
      sj=sjp;
      wj=wjp;
     end
    end
				
   elseif order==2
    if norm==1
     hj=h(kptk(ic)+1);
     for j=1:kncp(ic)
      if j<kncp(ic)-1
       hjp=h(kptk(ic)+j+1);
       wjp=wreg2(ic,j+1);
       vjp=(c(kcpt+j+2)-c(kcpt+j+1))/hjp-(c(kcpt+j+1)-c(kcpt+j))/hj;
       vjp=ocdabsm(vjp,rhoab,kabs);
       df(kcpt+j)=df(kcpt+j)+wjp*vjp/hj;
      end
      if (j<kncp(ic) & j>1)
       df(kcpt+j)=df(kcpt+j)-wj*(1/hjm+1/hj)*vj;
      end
      if j>2, df(kcpt+j)=df(kcpt+j)+wjm*vjm/hjm; end
      wjm=wj;
      wj=wjp;
      hjm=hj;
      hj=hjp;
      vjm=vj;
      vj=vjp;
     end
    elseif norm==2
     hj=h(kptk(ic)+1);
     for j=1:kncp(ic)
      if j<kncp(ic)-1
       hjp=h(kptk(ic)+j+1);
       wjp=wreg2(ic,j+1);
       vjp=(c(kcpt+j+2)-c(kcpt+j+1))/hjp-(c(kcpt+j+1)-c(kcpt+j))/hj;
       df(kcpt+j)=df(kcpt+j)+wjp*2*vjp/hj;
      end
      if (j<kncp(ic) & j>1), df(kcpt+j)=df(kcpt+j)-wj*(2/hjm+2/hj)*vj; end
      if j>2, df(kcpt+j)=df(kcpt+j)+wjm*2*vjm/hjm; end
      wjm=wj;
      wj=wjp;
      hjm=hj;
      hj=hjp;
      vjm=vj;
      vj=vjp;
     end
    end
   end		

  end
 end
end

% End of ocdregzn

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
