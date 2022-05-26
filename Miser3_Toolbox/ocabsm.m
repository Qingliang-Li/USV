function sabs = ocabsm(g,rhoab,kabs)
% The smoothing function for abs(g), namely abs(g) if abs(g)>=rhoab,
% else    abs(g) = (g*g+rho*rho)/(2*rho)  for  kabs=1 
%                = sqrt(g*g+rhoab*rhoab)  for  kabs=2 

if abs(g) >= rhoab
 sabs=abs(g);
elseif kabs==1
 sabs=(g*g+rhoab*rhoab)/(rhoab+rhoab);
else	
 sabs=sqrt(g*g+rhoab*rhoab);
end
