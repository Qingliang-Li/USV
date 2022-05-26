function dsabs = ocdabsm(g,rhoab,kabs)
% Derivative wrt argument of ocabsm.

if g >= rhoab
 dsabs=1;
elseif g<=-rhoab
 dsabs=-1;
elseif kabs==1
 dsabs=g/rhoab;
else	
 dsabs=g/sqrt(g*g+rhoab*rhoab);
end
