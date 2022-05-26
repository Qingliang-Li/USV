function srho = ocrhosm(g,rhoab)
% The rho smoothing of -min(0,g) as in Rehbock, Teo, Jennings, Lee 1989,
% for penalty formulations of all-time state constraints.                                          (g*g+rho*rho)/(2*rho)

if g >= rhoab
 srho=0;
elseif g<=-rhoab
 srho=-g;
else	
 srho=0.25*(g-rhoab)^2/rhoab;
end
