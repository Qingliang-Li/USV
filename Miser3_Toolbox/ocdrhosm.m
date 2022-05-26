function dsrho = ocdrhosm(g,rhoab)
% Derivative wrt argument of ocrhosm.

if g >= rhoab
 dsrho=0;
elseif g<=-rhoab
 dsrho=-1;
else	
 dsrho=0.5*(g-rhoab)/rhoab;
end
