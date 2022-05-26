function dy = ocstate(t,x,flag,zt,kint,kpar,koncts,tk,kptk,kpcp,upar)

global Cb

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

dy = ocf(t,x,ut,zt,upar);


