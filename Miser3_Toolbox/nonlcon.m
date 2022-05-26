function [cineq,ceq,grcineq,grceq]=nonlcon(c,ocsame,tolx,reltest,kpqk,tknot,...
		dtk,h,tk,tqd,tau,kpqtau,kpar,koncts,kpcp,kncp,kptk,ns,nz,lget,...
		numder,epsjt,taujt,rhoab,tolpsi,tfinal,ferr,kout,consdata1,misc,...
		kreg,wreg0,wreg1,wreg2,lreg);
% routine to compute constraints and constraint gradients for FMINCON

%	consdata1=[ng,ngueq,ngceq,ngzeq,nguin,ngcin,ngzin]

ngc=consdata1(3)+consdata1(6);
nct=length(c)-nz;
upar=misc(7:6+misc(6));
if nargout <= 2
 flag=0020;
 if ocsame, flag=0010; end

 [dum,dum,g]= ocfgeval(flag,c,tolx,reltest,kpqk,tknot,...
  dtk,tk,tqd,tau,kpqtau,kpar,koncts,kpcp,kptk,ns,nz,ngc,lget,numder,...
  epsjt,taujt,rhoab,tolpsi,tfinal,ferr,kout,misc);
		
 ceq=[];cineq=[];
		
% Rearrange inequality and equality constraints into order
% Move canonical equalities first.
 if consdata1(3) > 0
  ceq(1:consdata1(3))=g(1:consdata1(3));
 end
% Move canonical inequalities next. Change inequalities from >= to <=
 if consdata1(6) > 0
  cineq(1:consdata1(6))=-g(consdata1(3)+1:consdata1(3)+consdata1(6));
 end

% Now for system parameter only constraints
 if consdata1(4)+consdata1(7) > 0
  gz = ocgz(c(nct+1:end),upar);
% Move system parameter equalities next.
  if consdata1(4) > 0
   ceq(consdata1(3)+1:consdata1(3)+consdata1(4))=gz(1:consdata1(4));
  end
% Move system parameter inequalities next. Change inequalities from >= to <=
  if consdata1(7) > 0
   cineq(consdata1(6)+1:consdata1(6)+consdata1(7))=...
    -gz(consdata1(4)+1:consdata1(4)+consdata1(7));
  end
 end
		
else
 flag=2020;
 if ocsame, flag=1010; end
	
 [dum,dum,g,gradg]= ocfgeval(flag,c,tolx,reltest,kpqk,tknot,...
  dtk,tk,tqd,tau,kpqtau,kpar,koncts,kpcp,kptk,ns,nz,ngc,lget,numder,...
  epsjt,taujt,rhoab,tolpsi,tfinal,ferr,kout,misc);
		
 ceq=[];cineq=[];grceq=[];grcineq=[];
		
% Rearrange inequality and equality constraints into order
% Move canonical equalities first.
 if consdata1(3) > 0
  ceq(1:consdata1(3))=g(1:consdata1(3));
 end
% Move canonical inequalities next. Change inequalities from >= to <=
 if consdata1(6) > 0
  cineq(1:consdata1(6))=-g(consdata1(3)+1:consdata1(3)+consdata1(6));
 end

% Now for system parameter only constraints
 if consdata1(4)+consdata1(7) > 0
  gz = ocgz(c(nct+1:end),upar);
% Move system parameter equalities next.
  if consdata1(4) > 0
   ceq(consdata1(3)+1:consdata1(3)+consdata1(4))=gz(1:consdata1(4));
  end
% Move system parameter inequalities next. Change inequalities from >= to <=
  if consdata1(7) > 0
   cineq(consdata1(6)+1:consdata1(6)+consdata1(7))=...
    -gz(consdata1(4)+1:consdata1(4)+consdata1(7));
  end
 end
	
% Rearrange gradients of inequality and equality constraints into order

% Rearrange inequality and equality constraints into order
% Move canonical equalities first.
 if consdata1(3) > 0
  grceq(:,1:consdata1(3))=gradg(1:consdata1(3),:)';
 end
% Move canonical inequalities next. Change inequalities from >= to <=
 if consdata1(6) > 0
  grcineq(:,1:consdata1(6))=-gradg(consdata1(3)+1:consdata1(3)+consdata1(6),:)';
 end

% Now for system parameter only constraints
 if consdata1(4)+consdata1(7) > 0
  dgzgz = ocdgzdz(c(nct+1:end),upar);
% Move system parameter equalities first.
  if consdata1(4) > 0
   grceq(:,consdata1(3)+1:consdata1(3)+consdata1(4))=...
    [zeros(consdata1(4),nct) dgzdz(1:consdata1(4),:)]';
  end
% Move system parameter inequalities next. Change inequalities from >= to <=
  if consdata1(7) > 0
   grcineq(:,consdata1(6)+1:consdata1(6)+consdata1(7))=...
    [zeros(consdata1(7),nz) -dgzdz(consdata1(4)+1:consdata1(4)+consdata1(7))]';
  end
 end
	
end

% End of nonlcon

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
