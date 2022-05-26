function [tqd,xqd,uqd]=miser3

% Main function for Miser3 Matlab Version 3.2
% Last change: July 2002
clc
clear all 

[terminate,title,lcdate,typefi,version,rsinfo,sysdata,knotdata,...
 condata1,condata2,pardata,consdata1,consdata2,acctol,opt,numder,kout,...
 fileinout,misc] = mkdata;

if terminate
 tqd=[];xqd=[];uqd=[];

else

% unpack some of the problem variables
 ns=sysdata(1);nc=sysdata(2);nz=sysdata(3);
 tstart=sysdata(4);tfinal=sysdata(5);
 lab=sysdata(6);
 nksets=sysdata(7);
 ng=consdata1(1);ngueq=consdata1(2);ngceq=consdata1(3);
 nguin=consdata1(5);ngcin=consdata1(6);
 ngu=ngueq+nguin;ngc=ngceq+ngcin;
 tolx=acctol(1);tolpsi=acctol(2);
 hmax=acctol(3);
 reltest=acctol(4);
 epsjts=acctol(5);taujts=acctol(6);
 rhoabs=acctol(7);
 maxite=opt(1);maxfun=opt(2);
 optprt=opt(3);
 epsopt=opt(4);accuracy=epsopt;
 epscon=opt(5);
 epscons=epscon;
 ncheck=misc(1);
 kabs=misc(5);
 nupar=misc(6);
 upar=misc(7:6+nupar);
 i1=ng+ngu*(nc+1)+ngc;
 lget=consdata2(i1+1:i1+ngc);
 let=0;if sum(lget>0),let=1;end

% Call setup to set up all the data structures for this problem.
 [koncts,kpcp,kncp,cspar,clower,cupper,lreg,kreg,wreg0,wreg1,wreg2,h,...
  tk,dtk,tknot,kptk,kpqk,tqd,tau,kpqtau,kpar,Auin,Buin,Aueq,Bueq,ngui,ferr]=...
  setup(ns,nc,nz,tstart,tfinal,nksets,knotdata,condata1,condata2,pardata,...
  consdata1,consdata2,hmax,title,lcdate,typefi,version,rsinfo,fileinout);

% Set some global variables
 global Once Cb Logicvalues Ntimes
 Once=0;
 Cb=999*ones(size(cspar));
 Logicvalues=ones(1,7);
 Ntimes=0;
 global xqd dxqd uqd psiqd
 uqd=[];

% print the values of user's fixed parameters.
 if nupar > 0
  fprintf('\n\n User parameters for this run are:\n');
  fprintf(' %15.7e%15.7e%15.7e%15.7e%15.7e\n',upar);
  fprintf('\n\n');
 end
 tic  % Start the timing

% save the original control and system parameters in history.
 histcount=0;
 histpar(1,:)=cspar;
% set the epsilon-tau algorithm for inequality functional constraints.
 epsjt=epsjts;
 taujt=taujts;
% set up smoothing of absolute values of integrands and elsewhere.
 rhoab=rhoabs;

% call ocfgeval to obtain values of all quantities at the initial point
 [f0,gradf0,g,gradg]= ocfgeval(2222,cspar,tolx,reltest,kpqk,tknot,...
  dtk,tk,tqd,tau,kpqtau,kpar,koncts,kpcp,kptk,ns,nz,ngc,lget,...
  numder,epsjt,taujt,rhoab,tolpsi,tfinal,ferr,kout,misc);

%______________________________________________________________________________
% Call the optimization routine
% This is the start of the outer optimisation loop
% Note that, with fmincon, there is no difference between a cold and
% a warm start
 continueopt=1;
 while continueopt

  if let 
   fprintf('\n EPS and TAU for max(0,g) smoothing  %12.4e%12.4e \n',...
    epsjt,taujt);
   accuracy=max([epsopt,0.1*epsjt]);
   epscons=max([epscon,0.01*accuracy,0.01*taujt]);
  end
  if lab
   fprintf('\n RHO for absolute value smoothing %12.4e \n',rhoab);
  end

% Set up the call to FMINCON
  lsearch='quadcubic';
  if opt(6)==1, lsearch='cubicpoly'; end
  optout='off';
  if optprt==1
   optout='final';
  elseif optprt==2
   optout='iter';
  end
  ocsame=0;

  options=optimset('DerivativeCheck','off','Diagnostics','off',...
   'DiffMaxChange',1e-06,'DiffMinChange',1e-8,'Display', optout,...
   'GradConstr','on','GradObj','on','LargeScale','off','MaxFunEvals',...
   maxfun,'MaxIter',maxite-1,'TolCon',epscons,'TolFun',accuracy,...
   'TolX',accuracy);
% 'LineSearchType',lsearch);
		
  [cspar,g0,nfail,output]=fmincon('func',cspar,Auin,Buin,Aueq,Bueq,clower,...
   cupper,'nonlcon',options,ocsame,tolx,reltest,kpqk,tknot,dtk,...
   h,tk,tqd,tau,kpqtau,kpar,koncts,kpcp,kncp,kptk,ns,nz,lget,numder,...
   epsjt,taujt,rhoab,tolpsi,tfinal,ferr,kout,consdata1,misc,...
   kreg,wreg0,wreg1,wreg2,lreg);
		
  [conin,coneq]=nonlcon(cspar,ocsame,tolx,reltest,kpqk,tknot,dtk,...
   h,tk,tqd,tau,kpqtau,kpar,koncts,kpcp,kncp,kptk,ns,nz,lget,numder,...
   epsjt,taujt,rhoab,tolpsi,tfinal,ferr,kout,consdata1,misc,...
   kreg,wreg0,wreg1,wreg2,lreg);
  iterno=output.iterations;
  funcno=output.funcCount;

% nfail > 0		solution reached
% nfail = 0		max number of iterations or function evals reached
% nfail < 0		failed to converge

  if nguin > 0, conin=[(Buin-Auin*cspar')',conin]; end
  if ngueq > 0, coneq=[(Bueq-Aueq*cspar')',coneq]; end
  g=[conin,coneq];
	
%______________________________________________________________________________
% What to do when the optimization routine has done something.
% Take a history first, then save, then check sequential optimization.
% Note that the index of histpar is one higher than the other history variables
  histcount=histcount+1;
  histf(histcount)=g0;
  if consdata1(1) > 0
   histg(histcount,:)=g;
  else
   histg(histcount)=0;
  end
  histabs(histcount)=rhoab;
  histeps(histcount)=epsjt;
  histtau(histcount)=taujt;
  histpar(histcount+1,:)=cspar;
  nfailhist(histcount)=nfail;
  iterhist(histcount)=iterno;
% Save this solution
  if kout(4) > 0
   kount(histcount)=ocsave(1,histcount,epsjt,taujt,rhoab,lab,let,kpcp,...
    kncp,cspar,title,lcdate,typefi,version,rsinfo,sysdata,knotdata,...
	 condata1,condata2,pardata,consdata1,consdata2,acctol,opt,numder,...
	 kout,fileinout,misc);
  else
   kount(histcount)=0;
  end
				
% Sort out action to be taken on failure.
  choice=-1;
  continueopt=0;
  if nfail <= 0
% choice = 0, stop
%        = 1, no sequential opt - restart the opt
%        = 2, continue with sequential opt
   choice=failurechoice(nfail,let,lab); 
   if choice==1, continueopt=1; end
  end

%______________________________________________________________________________
% Sort out the sequential optimization progress
  if lab & (nfail > 0 | choice==2)
   if rhoab > 0.99e-3
    rhoab=rhoab*0.1;
    continueopt=1;
   elseif let==0
    continueopt=0;
   end
  end	
%______________________________________________________________________________
  if let & (choice==2 | nfail > 0)
   continueopt=1;
   lfeas=1;
   epstj=epsjt;
   epsjt=0.0;
   tautj=taujt;
   taujt=0.0;
% call ocfgeval with epsjt and taujt both zero
   [f,dum,gm,dum]= ocfgeval(0022,cspar,tolx,reltest,kpqk,tknot,dtk,...
    tk,tqd,tau,kpqtau,kpar,koncts,kpcp,kptk,ns,nz,ngc,lget,numder,...
    epsjt,taujt,rhoab,tolpsi,tfinal,ferr,kout,misc);
		
   epsjt=epstj;
   taujt=tautj;
   fprintf('\n');
   for ig=1:ngcin
    i=ngceq+ig;
    if lget(i)
     histg(histcount,i)=gm(i);f
     fprintf('\n Constraint using eps-tau is no.%3i. Value of Geps is%12.4e',...
      ig,g(i));
     fprintf('\n Value for feasibility test is%12.4e',gm(i));
     if gm(i) < -max([epscon,epsjt*0.0025*abs(tfinal-tstart)]), lfeas=0; end
    end
   end
   if lab, rhoab=epsjt;end
   if lfeas
    epsjt=epsjt*0.1;
    taujt=taujt*0.1;
   else
    taujt=taujt*sqrt(0.1);
   end
		
   if epsjt < 0.9e-4
    continueopt=0;
   elseif taujt < hmax*0.1*epsjt
    fprintf('\n It appears that at least one constraint using the');
    fprintf(' epsilon-tau algorithm');
    fprintf('\n cannot be made feasible. The canonical constraints are: ');
    for ig=1:ngcin
     i=ngceq+ig;
     if lget(i) & gm(i) < -max([epscon,epsjt*0.025*abs(tfinal-tstart)])
      fprintf('%3i',i);
     end
    end
    continueopt=0;
   end
	
   if continueopt & lfeas
    fprintf('\n Previous iterate feasible - decreased both eps and tau.');
   else
    fprintf('\n Previous iterate infeasible - decreased tau only.');
   end
  end
	
  if (let | lab) & continueopt
   fprintf('\n \n Next sequential optimization beginning:');			
  end

 end   % end of continueopt
	
%______________________________________________________________________________
% Save the final results with bells and whistles.
 fprintf('\n\n MISER3 HAS FINISHED THE OPTIMIZATION.\n');
 toc   % end of timing

 [f,gradf0,g,gradg]= ocfgeval(2222,cspar,tolx,reltest,kpqk,tknot,dtk,...
  tk,tqd,tau,kpqtau,kpar,koncts,kpcp,kptk,ns,nz,ngc,lget,numder,...
  epsjt,taujt,rhoab,tolpsi,tfinal,ferr,kout,misc);

 if kout(4) > 0
  kount(histcount)=ocsave(2,0,epsjt,taujt,rhoab,lab,let,kpcp,kncp,...
   cspar,title,lcdate,typefi,version,rsinfo,sysdata,knotdata,condata1,...
   condata2,pardata,consdata1,consdata2,acctol,opt,numder,kout,...
   fileinout,misc);
  fprintf('\n Finished writing the save file.\n');
 end
%______________________________________________________________________________
% Print the final results.
 gm=g;
 if let
  epstj=epsjt;
  epsjt=0.0;
  tautj=taujt;
  taujt=0.0;
%	call ocfgeval with epsjt and epstj both 0
  [f,gradf0,gm,gradg]= ocfgeval(0022,cspar,tolx,reltest,kpqk,tknot,dtk,...
   tk,tqd,tau,kpqtau,kpar,koncts,kpcp,kptk,ns,nz,ngc,lget,numder,...
   epsjt,taujt,rhoab,tolpsi,tfinal,ferr,kout,misc);
	
  epsjt=epstj;
  taujt=tautj;
 end

 ocprint(cspar,g0,f,g,gm,histcount,histf,histg,histabs,histeps,histtau,...
  histpar,nfailhist,iterhist,nfail,iterno,funcno,kount,Aueq,Bueq,...
  Auin,Buin,kncp,tknot,lget,kpqk,tqd,title,lcdate,version,...
  rsinfo,sysdata,condata1,consdata1,opt,kout,fileinout,misc,lreg);

%______________________________________________________________________________
% Restart file
 if kout(2) > 0
  kresta=0;
  if nc > 0, kresta=restartchoice; end
  ocrestart (kresta,kpcp,kncp,cspar,title,lcdate,typefi,rsinfo,...
   sysdata,knotdata,condata1,condata2,pardata,consdata1,consdata2,...
   acctol,opt,numder,kout,fileinout,misc)
 end

% Close the error file
 fclose(ferr);
		
%______________________________________________________________________________
% Plotting
 plotsoln(nc,tqd,tknot,kpqk,koncts);

end % of terminate

% End of miser3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [terminate,title,lcdate,typefi,version,rsinfo,sysdata,knotdata,...
 condata1,condata2,pardata,consdata1,consdata2,acctol,opt,numder,kout,...
 fileinout,misc] = mkdata

% Problem variables are stored in the following arrays
% system information:  
%   sysdata=[ns,nc,nz,tstart,tfinal,lab,nksets]
% knot set information:  
%   knotdata=[kntype(i),knset(i),atk(i,j)]; i=1:nksets j=1:1:knset(i)
% control information 1:  
%   condata1=[kkc(i),koncts(i),kreg(i,j=1:9)] as columns; i=1:nc
% control information 2:
%   condata2=[lwbnd(i,j),inival(i,j),upbnd(i,j),wreg0(i,j),wreg1(i,j),wreg2(i,j)];
%     i=1:nc j=1:kncp(i) )
% system parameter data:
%   pardata=[zlwbnd(i),zinival(i),zupbnd(i)] as columns; i=1,nz
% constraint data 1:
%   consdata1=[ng,ngueq,ngceq,ngzeq,nguin,ngcin,ngzin]
% constraint data 2:			ngu=ngueq+nguin  ngc=ngceq+ngcin
%   consdata2=[contype(i); i=1:ng alpha(1):alpha(nc) beta(1) ... 
%     alpha((ngu-1)*nc+1):alpha(ngu*nc) beta(ngu) tau(j);j=1:ngc ...
%     lget(j);j=1:ngc lgab(j);j=1:ngc]
% accuracy and tolerance:
%   acctol=[tolx,tolpsi,hmax,reltest,epsjts,taujts,rhoabs]
% optimization selection:
%   opt=[maxite,maxfun,optprt,epsopt,epscon,lsearch]
% user derivatives selection:
%   numder(i,j),i=1:3,j=1:5
% input, output, and error file info:
%   kout=kout(i); i=1:4
%   fileinout=[f1,f2,f3,f4]
% miscellaneous information:
%   misc=[ncheck ststiff costiff nsave kabs nupar upar(i); i=1:nupar]

% Set default figure window colour for Miser
set(0,'DefaultFigureColor',[0.2 0.5 1]);

% fname2 = name of new or existing problem data input file
% stat = 'new' or 'old': whether a new or old file
% ledit = 0 , 1: whether or not the old file is to be edited
[fname2,stat,ledit,terminate]=welcome;

if stat=='new'
 d=datenum(now);
 d1=datestr(d,'ddd');d2=datestr(d,'mmm');d3=datestr(d,'dd');
 d4=datestr(d,'HH:MM:SS');d5=datestr(d,'yyyy');
 lcdate=char(strcat(d1,{' '},d2,{' '},d3,{' '},d4,{' '},d5));
 editsection=ones(1,11);
 typefi='INPUT';
 version=' 2.0';
 rsinfo='EDIT0000, RESTART0000';
 title=[];sysdata=[];knotdata=[];condata1=[];condata2=[];pardata=[];
 consdata1=[];consdata2=[];acctol=[];opt=[];numder=[];kout=[];
 fileinout=[];misc=[];	
else
 [title,lcdate,typefi,version,rsinfo,sysdata,knotdata,condata1,...
  condata2,pardata,consdata1,consdata2,acctol,opt,numder,kout,...
  fileinout,misc] = readfile(fname2);
end
	
if ledit
 if stat=='old'
% Determine which sections of the existing file are to be edited
  hfig1 = figure('Name','Miser3.2: Edit Information',...
   'Menubar', 'none','NumberTitle','off'); 
  h11 = uicontrol('Style','text','FontUnits','points','FontSize',14, ...
   'Units','normalized','Position',[0.25 0.9 0.6 0.06],...
   'String','Indicate which sections you wish to edit');
  String1={'File title information','System information',...
   'Knot information','Control information',...
   'System parameter information','Constraint information',...
   'Accuracy information','Optimization information',...
   'User derivatives information','Input, output and error file names',...
   'Miscellaneous information (including user parameters)'};
  for i=1:11
   c(i) = uicontrol('Style','checkbox','Units','normalized',...
    'Position',[0.25 0.84-0.07*(i-1) 0.6 0.05],'String',String1{i},...
    'Value',0);
  end
  h12 = uicontrol('Style','push','Units','normalized',...
   'Position',[0.45 0.02 0.1 0.07],'Callback', 'uiresume(gcbf);',...
   'String','Done');
  uiwait(hfig1)
  for i=1:11
   editsection(i)= get(c(i), 'Value');
  end
  close(hfig1)
 end

% Create the new file or edit the existing one and return all the values
 [title,sysdata,knotdata,condata1,condata2,pardata,...
  consdata1,consdata2,acctol,opt,numder,kout,fileinout,misc]...
  = editdata(fname2,editsection,title,sysdata,knotdata,condata1,...
  condata2,pardata,consdata1,consdata2,acctol,opt,numder,kout,fileinout,misc);
	
 if stat=='old'
% Obtain the name of the edited data file
  hfig2 = figure('Name','Miser3.2: New (edited) file name',...
   'Menubar', 'none','NumberTitle','off');
  h21 = uicontrol('Style','text','FontUnits','points','FontSize',12, ...
   'Units','normalized','Position',[0.25 0.8 0.5 0.1],...
   'String','type in name of new (edited) file and then click ''Continue''');
  h22 = uicontrol('Style','edit','Units','normalized',...
   'Position',[0.35 0.7 0.3 0.06],'String',fname2);
  h23 = uicontrol('Style','push','Units','normalized',...
   'Position',[0.42 0.55 0.15 0.08],'Callback', 'uiresume(gcbf);',...
   'String','Continue');
  uiwait(hfig2)
  fname1=get(h22, 'String');
  d=datenum(now);
  d1=datestr(d,'ddd');d2=datestr(d,'mmm');d3=datestr(d,'dd');
  d4=datestr(d,'HH:MM:SS');d5=datestr(d,'yyyy');
  lcdate=char(strcat(d1,{' '},d2,{' '},d3,{' '},d4,{' '},d5));
  change=0;
  stat='unknown';
  jedit=sscanf(rsinfo(5:8),'%i');
  jrest=sscanf(rsinfo(18:21),'%i');
  jedit=jedit+1;
  rsinfo=sprintf('EDIT%4i, RESTART%4i                                       ',...
   jedit,jrest);
  if strcmp(fname1,fname2)==0
   message={'You have changed the data input file name';...
    'Do you want to change the output file names?'};
   h24 = uicontrol('Style','text','Units','normalized',...
    'FontUnits','points','FontSize',12,'Position',[0.2 0.35 0.6 0.1],...
    'String',message);
   h25 = uicontrol('Style','checkbox','Units','normalized',...
    'Position',[0.3 0.25 0.1 0.05],'String','Yes');
   h26 = uicontrol('Style','checkbox','Units','normalized',...
    'Position',[0.6 0.25 0.1 0.05],'String','No');
   h27 = uicontrol('Style','push','Units','normalized',...
    'Position',[0.4 0.15 0.12 0.07],'Callback','uiresume(gcbf);',...
    'String','Done');
   uiwait(hfig2);
   change= get(h25, 'Value');
   fname2=fname1;
  end	
  close(hfig2);
		
  if change==1
   editsection=[zeros(1,9) 1 0];
   [title,sysdata,knotdata,condata1,condata2,pardata,...
    consdata1,consdata2,acctol,opt,numder,kout,fileinout,misc]...
    = editdata(fname2,editsection,title,sysdata,knotdata,condata1,...
    condata2,pardata,consdata1,consdata2,acctol,opt,numder,kout,...
   fileinout,misc);
  end
 end

 fprintf('\n Writing the new data file \n');
 fidw=fopen(fname2,'w');
 writefile(fidw,title,lcdate,typefi,rsinfo,sysdata,knotdata,condata1,...
  condata2,pardata,consdata1,consdata2,acctol,opt,numder,kout,...
  fileinout,misc);
 fclose(fidw);

 if ~terminate
% Now read in edited file for problem solution
  [title,lcdate,typefi,version,rsinfo,sysdata,knotdata,condata1,...
   condata2,pardata,consdata1,consdata2,acctol,opt,numder,kout,...
   fileinout,misc] = readfile(fname2);
 end
		
end
	
% End of mkdata
		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [filename,stat,ledit,terminate]=welcome
% Function welcome obtains the name of an existing problem input data file
% or sets fnin2=' ', and determines if the file requires editing

hfig = figure('Name','Miser3.2: Optimal Control Software',...
 'Menubar', 'none','NumberTitle','off'); 

welcomemessage={'WELCOME TO MISER3.2 OPTIMAL CONTROL SOFTWARE';'  ';
 'Copyright held by'
 'L S Jennings and M E Fisher';'  ';
 'For information contact les@maths.uwa.edu.au'};
h11 = uicontrol('Style','text','FontUnits','points','FontSize',18,...
 'Units','normalized','Position',[0.1 0.83 0.8 0.13],'String',welcomemessage);

message2={' Choose one of the options below'}; 
h12 = uicontrol('Style','text','FontUnits','points','FontSize',14, ...
 'Units','normalized','Position',[0.2 0.72 0.6 0.05],'String',message2);
	
String={'Generate a new problem input data file',...
 'Read in and edit an existing data input file',...
 'Read in and execute an existing data input file'};

for i=1:3
 c(i) = uicontrol('Style','checkbox','Units','normalized',...
  'Position',[0.25 0.66-0.06*(i-1) 0.5 0.05],'String',String{i});
 end

h13 = uicontrol('Style','text','Units','normalized',...
 'Position',[0.15 0.43 0.7 0.05],'String',...
 'type in name of file (either new or existing) and then click ''Continue''' );		
h14 = uicontrol('Style','edit','Units','normalized',...
 'Position',[0.25 0.37 0.25 0.05]);	
h15 = uicontrol('Style','push','Units','normalized',...
 'Position',[0.65 0.37 0.15 0.05],'Callback', 'uiresume(gcbf);',...
 'String','Continue');

uiwait(hfig);

filename=get(h14, 'String');
for i=1:3
 choice(i)= get(c(i), 'Value');
end
ledit=1;
if choice(3) == 1, ledit=0; end
stat='old';
if choice(1) == 1, stat='new'; end

terminate=0;
if choice(3)==0
 message3={' Choose one of the options below and then click ''Done'''};
 h16 = uicontrol('Style','text','FontUnits','points','FontSize',14, ...
  'Units','normalized','Position',[0.1 0.26 0.8 0.05],'String',message3);
 String2={'Write the new data file and continue with MISER3',...
  'Write the new data file and terminate MISER3'};
 for i=1:2
  c2(i) = uicontrol('Style','checkbox','Units','normalized',...
   'Position',[0.25 0.2-0.06*(i-1) 0.6 0.05],'String',String2{i});
 end
 h17 = uicontrol('Style','push','Units','normalized',...
  'Position',[0.45 0.03 0.1 0.07],'Callback', 'uiresume(gcbf);',...
  'String','Done');	
 uiwait(hfig);
 terminate=get(c2(2), 'value');
end

close(hfig)
		
% End of welcome
		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [title,lcdate,typefi,version,rsinfo,sysdata,knotdata,condata1,...
 condata2,pardata,consdata1,consdata2,acctol,opt,numder,kout,...
 fileinout,misc] = readfile(filename)
	
% Reads in the problem input data file and places all the problem variables
% in temporary arrays

fid=fopen(filename,'r');

% Read in file information
fgetl(fid);
line=fgetl(fid);title=line(14:end);
line=fgetl(fid);lcdate=line(15:38);
line=fgetl(fid);typefi=line(13:end);
line=fgetl(fid);version=line(18:22);
line=fgetl(fid);rsinfo=line(16:end);

% Read in system information
fgetl(fid);fgetl(fid);
line=fgetl(fid);ns=sscanf(line(11:14),'%i');
line=fgetl(fid);nc=sscanf(line(11:14),'%i');
line=fgetl(fid);nz=sscanf(line(11:14),'%i');
line=fgetl(fid);tstart=sscanf(line(9:24),'%e');
line=fgetl(fid);tfinal=sscanf(line(9:24),'%e');
line=fgetl(fid);labch=sscanf(line(8),'%c');
lab=0;
if labch=='T', lab=1; end

% Read in knot sets definition
fgetl(fid);fgetl(fid);
line=fgetl(fid);nksets=sscanf(line(11:14),'%i');
sysdata=[ns,nc,nz,tstart,tfinal,lab,nksets];
fgetl(fid);
if nksets == 0
 kntype=[];
 knset=[];
 atk=[];
 fgetl(fid);fgetl(fid);
else	
 G1=fscanf(fid,'%i',[3 nksets]);fgetl(fid);
 G1a=G1';
 kntype=G1a(:,2);
 knset=G1a(:,3);
% Type 1 knot sets
 atk=zeros(nksets,max(knset));
 fgetl(fid);fgetl(fid);
 for i=1:nksets
  if kntype(i)==1
   fgetl(fid);
   G1b=fscanf(fid,'%e',[knset(i),1]);fgetl(fid);
   atk(i,1:knset(i))=G1b';
  else
   knotgap=(tfinal-tstart)/(knset(i)-1);
   atk(i,1:knset(i))=[tstart:knotgap:tfinal];
  end
 end
end		
knotdata=[kntype knset atk];

% Read in control definition
if nc > 0
 fgetl(fid);fgetl(fid);fgetl(fid);
 G2=fscanf(fid,'%i',[12 nc]);fgetl(fid);
 kkc=G2(2,:);
 koncts=G2(3,:);
 kreg=G2(4:12,:);
 condata1=[kkc' koncts' kreg'];

% Read in control lower, initial and upper values and regularization weights
 fgetl(fid);fgetl(fid);
 for ic=1:nc
  kncp(ic)=knset(kkc(ic))+koncts(ic)-1;
 end
 maxcol=max(kncp);
 condata2=zeros(nc,6*maxcol);
 for i=1:nc
  ncp=kncp(i);
  fgetl(fid);
  G3=fscanf(fid,'%i %e %e %e %e %e %e',[7 ncp]);fgetl(fid);
  condata2(i,1:ncp)=G3(2,:);
  condata2(i,maxcol+1:maxcol+ncp)=G3(3,:);
  condata2(i,2*maxcol+1:2*maxcol+ncp)=G3(4,:);
  condata2(i,3*maxcol+1:3*maxcol+ncp)=G3(5,:);
  condata2(i,4*maxcol+1:4*maxcol+ncp)=G3(6,:);
  condata2(i,5*maxcol+1:5*maxcol+ncp)=G3(7,:);
 end
	
else
 condata1=[];
 condata2=[];
end

% Read in system parameters
if nz > 0
 ncp=nc+1;
 fgetl(fid);fgetl(fid);fgetl(fid);
 G4=fscanf(fid,'%i %e %e %e',[4 nz]);fgetl(fid);
 zlwbnd=G4(2,:);
 zinival=G4(3,:);
 zupbnd=G4(4,:);
 pardata=[zlwbnd' zinival' zupbnd'];
else
 pardata=[];
end

% Read in constraint information
fgetl(fid);fgetl(fid);
G5a=fscanf(fid,'%i %i %i',[3 1]);fgetl(fid);
G5b=fscanf(fid,'%i %i %i',[3 1]);fgetl(fid);
consdata1(1)=sum(G5a)+sum(G5b);
consdata1(2:4)=G5a;
consdata1(5:7)=G5b;
ngu=consdata1(2)+consdata1(5);
ngc=consdata1(3)+consdata1(6);
ng=consdata1(1);
fgetl(fid);

if ng > 0
 consdata2=zeros(1,ng);
 for j=1:ng
  cline=fgetl(fid);
  contype(j)=sscanf(cline(7:8),'%i');
  taui(j)=sscanf(cline(23:37),'%e');
  set(j)=sscanf(cline(48),'%c');
  sab(j)=sscanf(cline(61),'%c');
  if contype(j)==1 | contype(j)==11
   G5c=fscanf(fid,'%e',[nc+1 1]);fgetl(fid);
   consdata2=[consdata2 G5c'];
  end
 end
 consdata2(1:ng)=contype;
 index=find(contype==2 | contype==12 | contype==13);
 tau=taui(index);
 tl=length(index);
 lget=zeros(1,tl);
 lgab=zeros(1,tl);
 for i=1:tl
  if set(index(i))=='y', lget(i)=1; end
  if sab(index(i))=='y', lgab(i)=1; end
 end
 consdata2=[consdata2 tau lget lgab];
   
else
 consdata2=[];
end

% Read in accuracy and tolerance
fgetl(fid);fgetl(fid);
line=fgetl(fid);tolx=sscanf(line(9:24),'%e');
line=fgetl(fid);tolpsi=sscanf(line(9:24),'%e');
line=fgetl(fid);hmax=sscanf(line(9:24),'%e');
line=fgetl(fid);reltest=sscanf(line(9:24),'%e');
line=fgetl(fid);epsjts=sscanf(line(9:24),'%e');
line=fgetl(fid);taujts=sscanf(line(9:24),'%e');
line=fgetl(fid);rhoabs=sscanf(line(9:24),'%e');
acctol=[tolx tolpsi hmax reltest epsjts taujts rhoabs];

% Read in optimization selection
fgetl(fid);fgetl(fid);
line=fgetl(fid);maxite=sscanf(line(11:14),'%i');
line=fgetl(fid);maxfun=sscanf(line(11:14),'%i');
line=fgetl(fid);optprt=sscanf(line(11:14),'%i');
line=fgetl(fid);epsopt=sscanf(line(9:24),'%e');
line=fgetl(fid);epscon=sscanf(line(9:24),'%e');
line=fgetl(fid);lsearch=sscanf(line(11:14),'%i');
opt=[maxite maxfun optprt epsopt epscon lsearch];

% Read in user derivatives selection
numder=zeros(3,5);
fgetl(fid);fgetl(fid);fgetl(fid);fgetl(fid);
for i=1:3
 line=fgetl(fid);
 lastj=4;
 if i==3, lastj=5; end
 for j=1:lastj
  if i~=2 | j~=3, numder(i,j)=sscanf(line(13+6*j),'%i'); end
 end
end

% Read in input, output and error file info
fgetl(fid);fgetl(fid);
line=fgetl(fid);kout(1)=sscanf(line(27:28),'%i');
f1=sscanf(line(39:end),'%c');
line=fgetl(fid);kout(2)=sscanf(line(27:28),'%i');
f2=sscanf(line(39:end),'%c');
line=fgetl(fid);kout(3)=sscanf(line(27:28),'%i');
f3=sscanf(line(39:end),'%c');
line=fgetl(fid);kout(4)=sscanf(line(27:28),'%i');
f4=sscanf(line(39:end),'%c');
fileinout=char(f1,f2,f3,f4);

% Read in miscellaneous selection
fgetl(fid);fgetl(fid);
line=fgetl(fid);ncheck=sscanf(line(11:14),'%i');
line=fgetl(fid);ststiff=sscanf(line(11:14),'%i');
line=fgetl(fid);costiff=sscanf(line(11:14),'%i');
line=fgetl(fid);nsave=sscanf(line(11:14),'%i');
line=fgetl(fid);kabs=sscanf(line(11:14),'%i');
line=fgetl(fid);nupar=sscanf(line(11:14),'%i');
misc=[ncheck ststiff costiff nsave kabs nupar];
if nupar > 0
 upar=fscanf(fid,'%e',[nupar 1]);fgetl(fid);
 misc=[misc upar'];
end

fclose(fid);
		
% End of readfile
		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [title,sysdata,knotdata,condata1,condata2,pardata,...
 consdata1,consdata2,acctol,opt,numder,kout,fileinout,misc]...
 = editdata(fname2,editsection,title,sysdata,knotdata,condata1,...
 condata2,pardata,consdata1,consdata2,acctol,opt,numder,kout,fileinout,misc);
	
% Obtains all the edited values from the user for the creation of a new
% problem input data file or edits an old problem input data file

Stringnum=cellstr(num2str([1:100]'));

%------------------------------------------------------------------------------
if editsection(1)~=0
% File information
 hfig0 = figure('Name','Miser3.2: File Information','Menubar','none',...
  'NumberTitle','off'); 
 h01 = uicontrol('Style','text','FontUnits','points','FontSize',12,...
  'Units','normalized','Position',[0.2 0.6 0.6 0.06], ...
  'String','Type in Problem Title and then click ''Done''' );
 h02 = uicontrol('Style','edit','Units','normalized',...
  'Position',[0.1 0.45 0.8 0.05]);
 h03 = uicontrol('Style','push','Units','normalized',...
  'Position',[0.44 0.3 0.12 0.07],'Callback','uiresume(gcbf);',...
  'String','Done');
 uiwait(hfig0);
 title=get(h02, 'String');
 close(hfig0);
end	

%------------------------------------------------------------------------------
if editsection(2)~=0
% system information: sysdata(1:6)
 hfig1 = figure('Name','Miser3.2: System Information','Menubar','none',...
  'NumberTitle','off'); 
 mess1a={' Edit, where necessary, the information below and'};
 mess1b={' then click ''Done'''};
 message1=strcat(mess1a,mess1b);
 
 finished=0;
 while ~finished
  h11 = uicontrol('Style','text','FontUnits','points','FontSize',14, ...
   'Units','normalized','Position',[0.2 0.7 0.6 0.12],'String',message1);	
  String1={'Number of state variables','Number of control variables',...
   'Number of system parameters','The initial time','The final time'};
  String2={'1';'1';'0';'0';'1'};
  for i=1:5
   c11(i) = uicontrol('Style','text','Units','normalized', ...
    'Position',[0.25 0.6-0.07*(i-1) 0.35 0.05],'String',String1{i});
   c12(i) = uicontrol('Style','edit','Units','normalized',...
    'Position',[0.65 0.6-0.07*(i-1) 0.1 0.05],'String',String2{i});
  end
  h12 = uicontrol('Style','checkbox','Units','normalized', ...
   'Position',[0.15 0.23 0.6 0.05],'String',...
   'Use absolute value smoothing for the objective function','Value',0);
  h13 = uicontrol('Style','push','Units','normalized',...
   'Position',[0.4 0.1 0.12 0.07],'Callback', 'uiresume(gcbf);',...
   'String','Done');
  uiwait(hfig1);

  for i=1:5
   sysdata(i)= str2num(get(c12(i), 'String'));
  end
  sysdata(6)=get(h12, 'Value');
  close(hfig1);
	
  if sysdata(1)<= 0 | sysdata(2)<0 | sysdata(3)<0
   hfig1 = figure('Name','Miser3.2: System Information','Menubar','none',...
	 'NumberTitle','off');
   message1a={'You must have a positive number of state variables and a';
    'nonnegative number of  control variables and system parameters'};
   h10 = uicontrol('Style','text','FontUnits','points','FontSize',12, ...
    'Units','normalized','Position',[0.05 0.88 0.9 0.1],'String',message1a);
   finished=0;
  else
   finished=1;
  end	
 end	
		
end

%------------------------------------------------------------------------------
if editsection(3)~=0
% knot set data: sysdata(7) and knotdata
 editsection(4)=1;	% edit controls if knot info changes
 hfig2 = figure('Name','Miser3.2: Knot Set Information','Menubar','none',...
  'NumberTitle','off'); 
	
 if sysdata(2)==0
% If no controls
  message2={'You have no control functions but you can define one knot set';
   'for piecewise integration of the differential equations';
   'Please check one box and then click ''Done''           '};
  h21 = uicontrol('Style','text','FontUnits','points','FontSize',10, ...
   'Units','normalized','Position',[0.1 0.84 0.8 0.12],'String',message2);	
  String3={'I want to define a knot set for integration of the ode',...
   'I do not want to define a knot set'}; 
  for i=1:2
   c21(i) = uicontrol('Style','checkbox','Units','normalized',...
    'Position',[0.1 0.75-0.07*(i-1) 0.55 0.05],'String',String3{i});
  end
  h22 = uicontrol('Style','push','Units','normalized',...
   'Position',[0.68 0.72 0.1 0.07],'Callback','uiresume(gcbf);',...
	'String','Done');			
  uiwait(hfig2)
	
  for i=1:2
   choice(i)= get(c21(i), 'Value');
  end
  sysdata(7)=0;
  if choice(1) == 1, sysdata(7)=1; end
  close(hfig2);
		
 elseif sysdata(2)==1
% If one control, then only one knot set required
  sysdata(7)=1;

 else
% Obtain number of knot sets if number of controls > 1
  h23 = uicontrol('Style','text','Units','normalized',...
   'Position',[0.1 0.7 0.35 0.07],...
   'String',' Specify number of different knot sets and then click ''Done''');
  h24 = uicontrol('Style','edit','Units','normalized',...
   'Position',[0.48 0.71 0.08 0.05],'String','1');
  h25 = uicontrol('Style','push','Units','normalized',...
   'Position',[0.6 0.7 0.1 0.07],'Callback', 'uiresume(gcbf);',...
   'String','Done');
  uiwait(hfig2);
  sysdata(7)=str2num(get(h24, 'String'));
  if sysdata(7) > sysdata(2), sysdata(7)=1; end		
 end

% Obtain the type of each knot set and numbers of knots in each knot set
 if sysdata(7)>0
  mess3a={'Knot Set             Number of Knots'};
  mess3b={'                               Knot Type'};
  message3=strcat(mess3a,mess3b);
  h26 = uicontrol('Style','text','FontUnits','points','FontSize',10, ...
   'Units','normalized','Position',[0.1 0.55 0.7 0.05],'String',message3);
  for i=1:sysdata(7)
   c22(i) = uicontrol('Style','text','Units','normalized', ...
    'Position',[0.14 0.48-(i-1)*0.06 0.06 0.05],...
    'String',Stringnum{i});
   c23(i) = uicontrol('Style','edit','Units','points','Units','normalized',...
    'Position',[0.33 0.48-(i-1)*0.06 0.08 0.05],'String','11');
   c24(i) = uicontrol('Style','checkbox','Units','normalized', ...
    'Position',[0.5 0.48-(i-1)*0.06 0.2 0.05],'String','equally spaced',...
	 'Value',1);
   c25(i) = uicontrol('Style','checkbox','Units','normalized', ...
    'Position',[0.71 0.48-(i-1)*0.06 0.2 0.05],'String','user specified',...
	 'Value',0);
  end
  h27 = uicontrol('Style','push','Units','normalized',...
   'Position',[0.4 0.1 0.1 0.07],'Callback', 'uiresume(gcbf);','String','Done');
  uiwait(hfig2);
		
  for i=1:sysdata(7)
   knset(i)=str2num(get(c23(i), 'String'));
   kntype(i)=get(c25(i), 'Value');
  end
  close(hfig2);
	
% Specify the knots for all knotsets
  atk=zeros(sysdata(7),max(knset));
  finished=0;
  while ~finished
   if any(kntype)
    hfig2a = figure('Name','Miser3.2: User Specified Knots',...
     'Menubar','none','NumberTitle','off');
    str1=sprintf('For each knot set, ''First Knot = tstart'' (%11.4e)',...
     sysdata(4));
    str2=sprintf(', the rest of knots are in increasing order');
    str3=sprintf(', ''Last Knot = tfinal'' (%11.4e)',sysdata(5));
    str4=strcat(str1,str2,str3);
    h281 = uicontrol('Style','text','FontUnits','points','FontSize',11,...
     'Units','normalized','Position',[0.05 0.89 0.9 0.1],'String',str4);	
   end
	
   ydist=0.9;
   for i=1:sysdata(7)
    if kntype(i)==1
     ydist=ydist-0.07;
     str5=sprintf('Specify the knots for knot set %2i (Total of %2i knots)',...
      i,knset(i));					
     h282 = uicontrol('Style','text','Units','normalized',...
      'Position',[0.2 ydist 0.6 0.05],'String',str5);
     numrows=ceil(knset(i)/8);
     str6=num2str(sysdata(4));
     for k=1:numrows
      imin=8*(k-1)+1;
      imax=min(8*k,knset(i));
      ydist=ydist-0.06;
      for j=imin:imax
       if j==1
        h283(j) = uicontrol('Style','edit','Units','normalized','Position',...
         [0.05+(j-imin)*0.12 ydist 0.1 0.05],'String',num2str(sysdata(4)));
       elseif j==knset(i)
        h283(j) = uicontrol('Style','edit','Units','normalized','Position',...
         [0.05+(j-imin)*0.12 ydist 0.1 0.05],'String',num2str(sysdata(5)));
       else
        h283(j) = uicontrol('Style','edit','Units','normalized','Position',...
         [0.05+(j-imin)*0.12 ydist 0.1 0.05]);
       end
      end
     end
     ydist=ydist-0.07;
     h284 = uicontrol('Style','push','Units','normalized',...
     'Position',[0.4 ydist 0.12 0.06],'Callback','uiresume(gcbf);',...
	  'String','Continue');
     uiwait(hfig2a);
     for j=1:knset(i)
      atk(i,j)=str2num(get(h283(j), 'String'));
     end						
    else
     knotgap=(sysdata(5)-sysdata(4))/(knset(i)-1);
     atk(i,1:knset(i))=[sysdata(4):knotgap:sysdata(5)];
    end
   end
   if any(kntype), close(hfig2a); end
   finished=1;
   for i=1:sysdata(7)
    for j=2:knset(i)
     if atk(i,j-1)>=atk(i,j), finished=0; end
    end
   end
   if ~finished
    hfig2b = figure('Name','Miser3.2: Error','Menubar','none',...
	  'NumberTitle','off');
    message2b={'There is an error in the knots you have provided.';
     'The knots are not in increasing order. Please enter them again.'};
    h10a = uicontrol('Style','text','FontUnits','points','FontSize',12, ...
     'Units','normalized','Position',[0.05 0.7 0.9 0.1],'String',message2b);
    h10b = uicontrol('Style','push','Units','normalized',...
     'Position',[0.4 0.5 0.1 0.07],'Callback', 'uiresume(gcbf);',...
     'String','OK');
    uiwait(hfig2b);
    close(hfig2b);
   end	
  end
	
  knotdata=[kntype' knset' atk];
 end
end

%------------------------------------------------------------------------------
if editsection(4)~=0 & sysdata(2) > 0
% control information: condata1

% lbtype(i) and ubtype(i) are set to 0
 hfig3 = figure('Name','Miser3.2: Control Information','Menubar','none',...
  'NumberTitle','off'); 
 h301 = uicontrol('Style','checkbox','Units','normalized', ...
  'Position',[0.1 0.93 0.6 0.06],'String',...
  '  Regularization of controls to be added to the objective','Value',0);
 h302 = uicontrol('Style','push','Units','normalized',...
  'Position',[0.75 0.93 0.12 0.06],'Callback', 'uiresume(gcbf);',...
  'String','Continue');
 uiwait(hfig3);
 lreg = get(h301, 'Value');
	
 message4a={'Control    Knot Set                   Control continuity'};
 message4b={'                      Type of control initial value'};
 message4=strcat(message4a,message4b);
 h31 = uicontrol('Style','text','Units','normalized',...
  'Position',[0.01 0.865 0.92 0.05],'String',message4);
 for i=1:sysdata(2)
  c31(i) = uicontrol('Style','text','Units','normalized', ...
   'Position',[0.01 0.81-(i-1)*0.055 0.06 0.05],'String',Stringnum{i});
  c32(i) = uicontrol('Style','edit','Units','normalized',...
   'Position',[0.1 0.81-(i-1)*0.055 0.06 0.05],'String','1');
  c33(i) = uicontrol('Style','checkbox','Units','normalized',...
   'Position',[0.18 0.81-(i-1)*0.055 0.24 0.05],...
   'String','piecewise constant','Value',1);
  c34(i) = uicontrol('Style','checkbox','Units','normalized',...
   'Position',[0.42 0.81-(i-1)*0.055 0.21 0.05],...
   'String','piecewise linear','Value',0);
  c35(i) = uicontrol('Style','checkbox','Units','normalized',...
   'Position',[0.63 0.81-(i-1)*0.055 0.13 0.05],'String','constant','Value',1);
  c36(i) = uicontrol('Style','checkbox','Units','normalized',...
   'Position',[0.76 0.81-(i-1)*0.055 0.24 0.05],...
   'String','piecewise constant','Value',0);
 end
 h32 = uicontrol('Style','push','Units','normalized',...
  'Position',[0.49 0.47 0.12 0.06],'Callback','uiresume(gcbf);',...
  'String','Continue');
 uiwait(hfig3);
   
 for i=1:sysdata(2)
  kkc(i)=str2num(get(c32(i), 'String'));	
  koncts(i)=get(c34(i), 'Value');
  ivtype(i)=get(c36(i), 'value');
  knset=knotdata(:,2);
  kncp(i)=knset(kkc(i))+koncts(i)-1;
 end
 condata1=[kkc' koncts'];

% control information: condata2
 lwbnd=zeros(sysdata(2),max(kncp));
 upbnd=zeros(sysdata(2),max(kncp));
 inival=zeros(sysdata(2),max(kncp));
 if ivtype == 0
  finished=0;
  while ~finished
   message5a={'        Control                Lower bound            '};
   message5b={'Upper bound                     Initial value     '};
   message5=strcat(message5a,message5b);
   h33 = uicontrol('Style','text','Units','normalized',...
    'Position',[0.02 0.395 0.9 0.05],'String',message5);
   for i=1:sysdata(2)
    c38(i) = uicontrol('Style','text','Units','normalized',...
     'Position',[0.1 0.34-(i-1)*0.055 0.06 0.05],'String',Stringnum{i});
    c39(i) = uicontrol('Style','edit','Units','normalized',...
     'Position',[0.22 0.34-(i-1)*0.055 0.2 0.05],'String','-1e+20');
    c391(i) = uicontrol('Style','edit','Units','normalized',...
     'Position',[0.45 0.34-(i-1)*0.055 0.2 0.05],'String',' 1e+20');
    c392(i) = uicontrol('Style','edit','Units','normalized',...
     'Position',[0.68 0.34-(i-1)*0.055 0.2 0.05],'String','0');
   end
   h34 = uicontrol('Style','push','Units','normalized',...
    'Position',[0.5 0.00 0.1 0.06],'Callback', 'uiresume(gcbf);',...
	 'String','Done');
   uiwait(hfig3);
		
   finished=1;
   for i=1:sysdata(2)
    lwbnd(i,1:kncp(i))=str2num(get(c39(i),'String'))*ones(1,kncp(i));	
    upbnd(i,1:kncp(i))=str2num(get(c391(i),'String'))*ones(1,kncp(i));
    inival(i,1:kncp(i))=str2num(get(c392(i),'String'))*ones(1,kncp(i));
    if any(lwbnd(i,:)>upbnd(i,:)), finished=0; end
    if any(inival(i,:)>upbnd(i,:)), finished=0; end
    if any(inival(i,:)<lwbnd(i,:)), finished=0; end
   end	
   close(hfig3);
   if ~finished
    hfig3 = figure('Name','Miser3.2: Control Information',...
     'Menubar','none','NumberTitle','off'); 
    message5c={'The upper bound must be less than the lower bound and';
     'the initial value must be between the lower and upper bounds';
     'Please enter these again'};
    h34a = uicontrol('Style','text','Units','normalized',...
     'Position',[0.05 0.7 0.9 0.12],'String',message5c);
   end
  end
			
 else
  close(hfig3);
  message6={'Knot Number             Control Initial Value'};
  for i=1:sysdata(2)
   finished=0;
   while ~finished
    hfig4(i) = figure('Name','Miser3.2: Control Values',...
     'Menubar','none','NumberTitle','off');
    h41(i) = uicontrol('Style','text','FontUnits','points',...
     'FontSize',12,'Units','normalized','Position',[0.3 0.945 0.28 0.05],...
     'String','Control Number');
    h42(i) = uicontrol('Style','text','FontUnits','points',...
     'FontSize',12,'Units','normalized','Position',[0.6 0.945 0.06 0.05],...
     'String',Stringnum{i});
    c41(i) = uicontrol('Style','text','Units','normalized',...
     'Position',[0.06 0.885 0.18 0.05],'String','Lower bound =');
    c42(i) = uicontrol('Style','edit','Units','normalized',...
     'Position',[0.26 0.885 0.2 0.05],'String','-1e+20');
    c43(i) = uicontrol('Style','text','Units','normalized',...
     'Position',[0.52 0.885 0.18 0.05],'String','Upper bound =');
    c44(i) = uicontrol('Style','edit','Units','normalized',...
     'Position',[0.72 0.885 0.2 0.05],'String',' 1e+20');
    if kncp(i) < 15
     h43(i) = uicontrol('Style','text','Units','normalized',...
       'Position',[0.27 0.825 0.42 0.05],'String',message6);
     for j=1:kncp(i)
      c45(j) = uicontrol('Style','text','Units','normalized',...
       'Position',[0.32 0.77-(j-1)*0.051 0.06 0.05],'String',Stringnum{j});
      c46(j) = uicontrol('Style','edit','Units','normalized',...
       'Position',[0.48 0.77-(j-1)*0.051 0.2 0.05],'String','0');
     end
    else
     h43(i) = uicontrol('Style','text','Units','normalized',...
      'Position',[0.05 0.825 0.42 0.05],'String',message6);
     h44(i)= uicontrol('Style','text','Units','normalized',...
      'Position',[0.53 0.825 0.42 0.05],'String',message6);
     for j=1:14
      c45(j) = uicontrol('Style','text','Units','normalized',...
       'Position',[0.08 0.77-(j-1)*0.051 0.06 0.05],'String',Stringnum{j});
      c46(j) = uicontrol('Style','edit','Units','normalized',...
       'Position',[0.26 0.77-(j-1)*0.051 0.2 0.05],'String','0');
     end
     for j=15:kncp(i)
      c45(j) = uicontrol('Style','text','Units','normalized',...
       'Position',[0.56 0.77-(j-15)*0.051 0.06 0.05],'String',Stringnum{j});
      c46(j) = uicontrol('Style','edit','Units','normalized',...
       'Position',[0.74 0.77-(j-15)*0.051 0.2 0.05],'String','0');
     end
    end
			
    h45(i) = uicontrol('Style','push','Units','points',...
     'Units','normalized','Position',[0.4 0.01 0.1 0.07],...
     'Callback', 'uiresume(gcbf);','String','Done');
    uiwait(hfig4(i));
    lwbnd(i,1:kncp(i))=str2num(get(c42(i), 'String'))*ones(1,kncp(i));	
    upbnd(i,1:kncp(i))=str2num(get(c44(i), 'String'))*ones(1,kncp(i));
    for j=1:kncp(i)
     inival(i,j)=str2num(get(c46(j), 'String'));
    end	
    close(hfig4(i));
    finished=1;
    if any(lwbnd(i,:)>upbnd(i,:)), finished=0; end
    if any(inival(i,:)>upbnd(i,:)), finished=0; end
    if any(inival(i,:)<lwbnd(i,:)), finished=0; end
    if ~finished
     hfig4a = figure('Name','Miser3.2: Error','Menubar','none',...
      'NumberTitle','off');
     message4d={'There is an error in the values you have provided.';
      'The upper bound must be less than the lower bound and';
      'the initial value must be between the lower and upper bounds.';
      'Please enter these again.'};
     h10a = uicontrol('Style','text','FontUnits','points','FontSize',12, ...
      'Units','normalized','Position',[0.05 0.7 0.9 0.2],'String',message4d);
     h10b = uicontrol('Style','push','Units','normalized',...
      'Position',[0.4 0.5 0.1 0.07],'Callback', 'uiresume(gcbf);',...
      'String','OK');
     uiwait(hfig4a);
     close(hfig4a);
    end
   end
  end
 end
	
 condata2=[lwbnd inival upbnd];
	
% Regularization: append kreg to condata1 and wreg0/1/2 to condata2

% First set default values
 kreg=-ones(sysdata(2),9);
 wreg0=zeros(sysdata(2),max(kncp));
 wreg1=zeros(sysdata(2),max(kncp));
 wreg2=zeros(sysdata(2),max(kncp));
		
 if lreg
  hfig46 = figure('Name','Miser3.2: Regularization Information',...
   'Menubar','none','NumberTitle','off'); 
  h464 = uicontrol('Style','text','Units','normalized', ...
   'Position',[0.05 0.84 0.65 0.05],'String',...
   'Do you want all control functions equally in the regularization?');
  h465 = uicontrol('Style','checkbox','Units','normalized',...
   'Position',[0.75 0.84 0.1 0.05],'String','Yes','Value',1);
  h466 = uicontrol('Style','checkbox','Units','normalized',...
   'Position',[0.88 0.84 0.1 0.05],'String','No','Value',0);
  h467 = uicontrol('Style','push','Units','normalized',...
   'Position',[0.4 0.77 0.12 0.06],'Callback', 'uiresume(gcbf);',...
   'String','Continue');
  uiwait(hfig46);
  default = get(h465, 'Value');
		
  if default
   string4a={'order (0, 1 or 2)      norm (1 or 2)      weight (0, 1 or 2)'};
   string4b={'0','1','0'};
   h468 = uicontrol('Style','text','Units','normalized',...
    'Position',[0.2 0.55 0.6 0.06],'String',...
    'Specify 3 digits for type of regularization');
   c461 = uicontrol('Style','text','Units','normalized', ...
    'Position',[0.1 0.48 0.8 0.05],'String',string4a);
   for i=1:3
    c462(i) = uicontrol('Style','edit','Units','normalized',...
     'Position',[0.25+(i-1)*0.2 0.41 0.1 0.05],'String',string4b{i});
   end
   h469 = uicontrol('Style','text','Units','normalized',...
    'Position',[0.2 0.26 0.4 0.05],'String','Specify common penalty weight');
   c463 = uicontrol('Style','edit','Units','normalized',...
    'Position',[0.61 0.26 0.1 0.05],'String','0');						
   h4691 = uicontrol('Style','push','Units','normalized',...
    'Position',[0.4 0.1 0.12 0.07],'Callback', 'uiresume(gcbf);',...
    'String','Done');
   uiwait(hfig46);
   for i=1:3
    regtype(i) = str2num(get(c462(i), 'String'));
   end
   defwt= str2num(get(c463, 'String'));
   kreg(:,regtype(1)*3+1:regtype(1)*3+3)=ones(sysdata(2),1)*regtype;
   if regtype(1)==0, wreg0=defwt*ones(sysdata(2),max(kncp)); end
   if regtype(1)==1, wreg1=defwt*ones(sysdata(2),max(kncp)); end
   if regtype(1)==2, wreg2=defwt*ones(sysdata(2),max(kncp)); end
			
   close(hfig46);
  else
   close(hfig46);
   for ic=1:sysdata(2)				
    hfig47(ic) = figure('Name','Miser3.2: Regularization Information',...
     'Menubar','none','NumberTitle','off');
    str1=sprintf('Regularization of control %2i added to objective?',ic);
    h471 = uicontrol('Style','text','FontUnits','points',...
     'FontSize',12,'Units','normalized',...
     'Position',[0.05 0.91 0.6 0.05],'String',str1);
    c471=uicontrol('Style','checkbox','Units','normalized', ...
     'Position',[0.68 0.91 0.08 0.05],'String','No');
    c472=uicontrol('Style','checkbox','Units','normalized', ...
     'Position',[0.78 0.91 0.09 0.05],'String','Yes','Value',1);
    h472 = uicontrol('Style','push','Units','points',...
     'Units','normalized','Position',[0.89 0.91 0.1 0.05],...
     'Callback', 'uiresume(gcbf);','String','Continue');
    uiwait(hfig47(ic));
    addcontrol=get(c472, 'Value');
    if addcontrol
     string4a={'order (0, 1 or 2)      norm (1 or 2)      weight (0, 1 or 2)'};
     string4b={'0','1','0'};
     h473 = uicontrol('Style','text','Units','normalized',...
      'Position',[0.2 0.8 0.6 0.05],'String',...
      'Specify 3 digits for type of regularization');
     c473 = uicontrol('Style','text','Units','normalized', ...
      'Position',[0.1 0.73 0.8 0.05],'String',string4a);
     for j=1:3
      c474(j) = uicontrol('Style','edit','Units','normalized',...
       'Position',[0.25+(j-1)*0.2 0.67 0.1 0.05],'String',string4b{j});
     end
     h474 = uicontrol('Style','push','Units','normalized',...
      'Position',[0.4 0.57 0.12 0.07],'Callback', 'uiresume(gcbf);',...
      'String','Continue');
     uiwait(hfig47(ic));
     for j=1:3
      regtype(j) = str2num(get(c474(j), 'String'));
     end
     kreg(ic,regtype(1)*3+1:regtype(1)*3+3)=regtype;
		
     if regtype(3)~=1
      str2=sprintf('Specify common penalty weight for regularization of order %2i',...
       regtype(1));
      h475 = uicontrol('Style','text','Units','normalized',...
       'Position',[0.1 0.46 0.65 0.05],'String',str2);
      c475 = uicontrol('Style','edit','Units','normalized',...
       'Position',[0.8 0.46 0.1 0.05],'String','0');						
      h476 = uicontrol('Style','push','Units','normalized',...
       'Position',[0.4 0.3 0.12 0.07],'Callback','uiresume(gcbf);',...
       'String','Done');
      uiwait(hfig47(ic));	
      defwt= str2num(get(c475, 'String'));
      if regtype(1)==0, wreg0(ic,:)=defwt*ones(1,max(kncp)); end
      if regtype(1)==1, wreg1(ic,:)=defwt*ones(1,max(kncp)); end
      if regtype(1)==2, wreg2(ic,:)=defwt*ones(1,max(kncp)); end
     else
      if regtype(1)==0
       i1=1; i2=kncp(ic);
      elseif regtype(1)==1 
       i1=2; i2=kncp(ic);
      else
       i1=2; i2=kncp(ic)-1;
      end
      str3=sprintf('Give list of %2i penalty weights for regularization of order %2i',...
       i2-i1+1,regtype(1));
      h477 = uicontrol('Style','text','Units','normalized',...
       'Position',[0.1 0.45 0.8 0.05],'String',str3);
      numrows=ceil((i2-i1+1)/8);
      for k=1:numrows
       imin=8*(k-1)+1;
       imax=min(8*k,i2-i1+1);
       for j=imin:imax
        h478(j) = uicontrol('Style','edit','Units','normalized',...
         'Position',[0.05+(j-imin)*0.12 0.39-(k-1)*0.07 0.1 0.05]);
       end
      end
      h479 = uicontrol('Style','push','Units','normalized',...
       'Position',[0.4 0.05 0.12 0.07],'Callback','uiresume(gcbf);',...
       'String','Done');
      uiwait(hfig47(ic));
      for j=i1:i2
       wt(j)=str2num(get(h478(j-i1+1), 'String'));
      end
      if regtype(1)==0, wreg0(ic,i1:i2)=wt(i1:i2); end
      if regtype(1)==1, wreg1(ic,i1:i2)=wt(i1:i2); end
      if regtype(1)==2, wreg2(ic,i1:i2)=wt(i1:i2); end
     end
    end
    close(hfig47(ic));
   end
  end	

% Set the regularization weights for all cases
  atk=knotdata(:,3:end);
  kptkset(1)=1;
  for ik=2:sysdata(7)
   kptkset(ik)=kptkset(ik-1)+knset(ik-1);
  end
  kptk=kptkset(kkc);
  tk=[];h=[];
  for i=1:sysdata(7)
   h(kptkset(i))=-1e20;
   h(kptkset(i)+1:kptkset(i)+knset(i)-1)=atk(i,2:knset(i))-atk(i,1:knset(i)-1);
  end
  for ic=1:sysdata(2)
   regtype=[-1 -1 -1];			
   [j1 j2 findreg]=find(kreg(ic,:)+1);
   if length(findreg)>0, regtype=findreg-1; end
   kp=kptk(ic);kn=kncp(ic);
   if regtype(1)==0 & regtype(3)==2
    if koncts(ic)==0
     wreg0(ic,1:kn)=h(kp+1:kp+kn).*wreg0(ic,1:kn);
    elseif koncts(ic)==1
     wreg0(ic,1)=1e20;
     hi=h(kp+1:kp+kn-1);
     if regtype(2)==1, hi=hi/2; end
     if regtype(2)==2, hi=hi/3; end
     wreg0(ic,2:kn)=hi.*wreg0(ic,2:kn);
    end
   elseif regtype(1)==1
    wreg1(ic,1)=1e20;
    if regtype(2)==2  & regtype(3)==2
     wreg1(ic,2:kn)=wreg1(ic,2:kn)./h(kp+1:kp+kn-1);
    end	
   elseif regtype(1)==2
    wreg2(ic,1)=1e20;
    wreg2(ic,kn)=1e20;
    hi=h(kp+1:kp+kn-2);
    hip=h(kp+2:kp+kn-1);
    if (regtype(2)==1 & regtype(3)==1) | (regtype(2)==2 & regtype(3)==2)
     wreg2(ic,2:kn-1)=2*wreg2(ic,2:kn-1)./(hi+hip);
    end
    if (regtype(2)==2 & regtype(3) < 2)
     wreg2(ic,2:kn-1)=4*wreg2(ic,2:kn-1)./(hi+hip).^2;
    end
   end
  end
		
 end
	
 condata1=[condata1 kreg];
 condata2=[condata2 wreg0 wreg1 wreg2];
	
end

%------------------------------------------------------------------------------
if editsection(5)~=0 & sysdata(3) > 0
% system parameter information: pardata

 finished=0;
 while ~finished
  hfig5 = figure('Name','Miser3.2: System Parameter Information',...
   'Menubar','none','NumberTitle','off'); 
  message7a={'System Parameter           Lower Bound          '};
  message7b={'Initial Value          Upper Bound'};
  message7=strcat(message7a,message7b);

  h51 = uicontrol('Style','text','Units','normalized',...
   'Position',[0.02 0.88 0.76 0.05],'String',message7);
  for i=1:sysdata(3)
   c51(i) = uicontrol('Style','text','Units','normalized',...
    'Position',[0.1 0.8-(i-1)*0.055 0.06 0.05],'String',Stringnum{i});
   c52(i) = uicontrol('Style','edit','Units','normalized',...
    'Position',[0.3 0.8-(i-1)*0.055 0.1 0.05],'String','-100');
   c53(i) = uicontrol('Style','edit','Units','normalized',...
    'Position',[0.48 0.8-(i-1)*0.055 0.1 0.05],'String','0');
   c54(i) = uicontrol('Style','edit','Units','normalized',...
    'Position',[0.66 0.8-(i-1)*0.055 0.1 0.05],'String','100');
  end
		
  h52 = uicontrol('Style','push','Units','normalized',...
   'Position',[0.4 0.2 0.1 0.07],'Callback', 'uiresume(gcbf);',...
   'String','Done');
  uiwait(hfig5);
  for i=1:sysdata(3)
   zlwbnd(i)=str2num(get(c52(i), 'String'));	
   zinival(i)=str2num(get(c53(i), 'String'));
   zupbnd(i)=str2num(get(c54(i), 'String'));
  end
  close(hfig5);
	
  finished=1;
  if any(zlwbnd>zupbnd), finished=0; end
  if any(zinival>zupbnd), finished=0; end
  if any(zinival<zlwbnd), finished=0; end
  if ~finished
   hfig5a = figure('Name','Miser3.2: System Parameter Information',...
    'Menubar','none','NumberTitle','off');
   message5d={'There is an error in the values you have provided.';
    'The upper bound must be less than the lower bound and';
    'the initial value must be between the lower and upper bounds.';
    'Please enter these again.'};
   h51a = uicontrol('Style','text','FontUnits','points','FontSize',12, ...
    'Units','normalized','Position',[0.05 0.7 0.9 0.2],'String',message5d);
   h51b = uicontrol('Style','push','Units','normalized',...
    'Position',[0.4 0.5 0.1 0.07],'Callback', 'uiresume(gcbf);',...
    'String','OK');
   uiwait(hfig5a);
   close(hfig5a);
  end
 end % finished

 pardata=[zlwbnd' zinival' zupbnd'];

end

%------------------------------------------------------------------------------
if editsection(6)~=0
% constraint information - numbers of constraints: consdata1

 hfig6 = figure('Name','Miser3.2: Numbers of Constraints',...
  'Menubar','none','NumberTitle','off'); 
	
 Stringc1={'Number of linear control equality constraints',...
  'Number of canonical equality constraints',...
  'Number of system parameter equality constraints',...
  'Number of linear control inequality constraints',...
  'Number of canonical inequality constraints',...
  'Number of system parameter inequality constraints'};
 for i=1:6
  c61(i) = uicontrol('Style','text','Units','normalized', ...
   'Position',[0.1 0.75-0.08*(i-1) 0.55 0.05],'String',Stringc1{i});
  c62(i) = uicontrol('Style','edit','Units','normalized',...
   'Position',[0.7 0.75-0.08*(i-1) 0.1 0.05],'String','0');
 end
 h61 = uicontrol('Style','push','Units','normalized',...
  'Position',[0.4 0.2 0.12 0.07],'Callback', 'uiresume(gcbf);','String','Done');
 uiwait(hfig6);

 for i=1:6
  consdata1(i+1)= str2num(get(c62(i), 'String'));
 end
 consdata1(1)=sum(consdata1(2:7));

 close(hfig6);
	
 contype=[];
 for i=1:consdata1(1)
  if i<=consdata1(2)
   contype(i)=1;
  elseif i<=sum(consdata1(2:3))
   contype(i)=2;
  elseif i<=sum(consdata1(2:4))
   contype(i)=3;
  elseif i<=sum(consdata1(2:5))
   contype(i)=11;
  elseif i<=sum(consdata1(2:6))
   contype(i)=12;
  else
   contype(i)=14;
  end
 end
   
% constraint information - constraint properties: consdata2
 ngu=consdata1(2)+consdata1(5);
 ngc=consdata1(3)+consdata1(6);
 alphabeta=[];tau=[];lgab=[];lget=[];
 contype=[];
   
 if ngu+ngc > 0   
  hfig7 = figure('Name','Miser3.2: Constraint Properties',...
   'Menubar','none','NumberTitle','off');

  if ngu > 0
   message8a=...
    {'Linear control constraints are in matrix form ''Au+b>=0'' with equality';
    'constraints first. Edit entries of the matrix A and the column vector b'};
   h71 = uicontrol('Style','text','Units','normalized',...
    'Position',[0.05 0.9 0.9 0.08],'String',message8a);
   h72 = uicontrol('Style','text','Units','normalized',...
    'Position',[0.05 0.8 0.08 0.05],'String','A = ');
   h73 = uicontrol('Style','text','Units','normalized',...
    'Position',[0.77 0.8 0.08 0.05],'String','b = ');
   for i=1:ngu
    for j=1:sysdata(2)
     h74(i,j) = uicontrol('Style','edit','Units','normalized',...
      'Position',[0.14+(j-1)*0.1 0.83-(i-1)*0.05 0.09 0.05],'String','0');
    end     
    h75(i) = uicontrol('Style','edit','Units','normalized',...
     'Position',[0.86 0.83-(i-1)*0.05 0.1 0.05],'String','0');
   end
   h76 = uicontrol('Style','push','Units','normalized',...
    'Position',[0.45 0.5 0.12 0.07],'Callback', 'uiresume(gcbf);',...
    'String','Done');
   uiwait(hfig7);
   for i=1:ngu
    for j=1:sysdata(2)
     alpha(j)=str2num(get(h74(i,j), 'String'));
    end         
    beta=str2num(get(h75(i), 'String'));
    alphabeta=[alphabeta alpha beta];
   end
  end
   
  if ngc > 0
   mess8b1={'Canonical constraints: equality constraints'};
   mess8b2={' before inequality constraints'};
   message8b=strcat(mess8b1,mess8b2);
   h77 = uicontrol('Style','text','Units','normalized',...
    'Position',[0.1 0.41 0.8 0.05],'String',message8b);
   for i=1:ngc
    c71(i) = uicontrol('Style','text','Units','normalized',...
     'Position',[0.01 0.34-(i-1)*0.05 0.05 0.05],'String',Stringnum{i});
    c72(i) = uicontrol('Style','text','Units','normalized',...
     'Position',[0.07 0.34-(i-1)*0.05 0.2 0.05],...
     'String','characteristic time:');
    c73(i) = uicontrol('Style','edit','Units','normalized',...
     'Position',[0.27 0.34-(i-1)*0.05 0.11 0.05],'String',num2str(sysdata(5)));
    c74(i) = uicontrol('Style','checkbox','Units','normalized',...
     'Position',[0.4 0.34-(i-1)*0.05 0.3 0.05],'String',...
     'absolute value smoothing');
    if i>consdata1(3)
     c75(i) = uicontrol('Style','checkbox','Units','normalized',...
      'Position',[0.72 0.34-(i-1)*0.05 0.27 0.05],'String',...
      'epsilon-tau smoothing');
    end
   end
   h88 = uicontrol('Style','push','Units','normalized',...
    'Position',[0.45 0.01 0.12 0.07],'Callback', 'uiresume(gcbf);',...
    'String','Done');
   uiwait(hfig7);

   for i=1:ngc
    tau(i)=str2num(get(c73(i), 'String'));      
    lgab(i)=get(c74(i), 'Value');
    lget(i)=0;
    if i>consdata1(3)
     lget(i)=get(c75(i), 'Value');
    end
   end
   if lgab~=0,sysdata(6)=1;end	%set lab='true'
			
  end
  
  close(hfig7);
	
  ig=0;
  for i=1:consdata1(1)
   if i<=consdata1(2)
    contype(i)=1;
   elseif i<=sum(consdata1(2:3))
    contype(i)=2;
    ig=ig+1;
   elseif i<=sum(consdata1(2:4))
    contype(i)=3;
   elseif i<=sum(consdata1(2:5))
    contype(i)=11;
   elseif i<=sum(consdata1(2:6))
    contype(i)=12;
    ig=ig+1;
    if lget(ig), contype(i)=13; end
   else
    contype(i)=14;
   end
  end
		
 end
	
 consdata2=[contype alphabeta tau lget lgab];

end

%------------------------------------------------------------------------------
if editsection(7)~=0
% accuracy and tolerance information: acctol

 hfig8 = figure('Name','Miser3.2: Accuracy and Tolerance',...
  'Menubar','none','NumberTitle','off');
 message81={' Edit, where necessary, the information below';
  'and then click ''Continue''';
  'The values shown are the default values'}; 
 h81 = uicontrol('Style','text','FontUnits','points','FontSize',12, ...
  'Units','normalized','Position',[0.2 0.78 0.6 0.16],'String',message81);	
 String81={'State equation tolerance','Costate equation tolerance',...
  'Quadrature interval relative maximum','Relative error test',...
  'Initial eps for eps-tau smoothing',...
  'Initial rho for abs smoothing'};
 String82={'1e-09';'1e-09';'1e-02';'1e+00';'1e-02';'1e-02'};
 for i=1:6
  c81(i) = uicontrol('Style','text','Units','normalized', ...
   'Position',[0.2 0.7-0.07*(i-1) 0.4 0.05],'String',String81{i});
  c82(i) = uicontrol('Style','edit','Units','normalized',...
   'Position',[0.65 0.7-0.07*(i-1) 0.15 0.05],'String',String82{i});
 end
 h82 = uicontrol('Style','push','Units','normalized',...
  'Position',[0.4 0.25 0.12 0.07],'Callback', 'uiresume(gcbf);',...
  'String','Continue');
 uiwait(hfig8);

 for i=1:5
  acctol(i)= str2num(get(c82(i), 'String'));
 end
 acctol(7) = str2num(get(c82(6), 'String'));
	
 taujts=acctol(5)*abs(sysdata(5)-sysdata(4))/4;
 c83 = uicontrol('Style','text','Units','normalized', ...
  'Position',[0.2 0.12 0.4 0.05],'String','Initial tau for eps-tau smoothing');
 c84 = uicontrol('Style','edit','Units','normalized',...
  'Position',[0.65 0.12 0.15 0.05],'String',num2str(taujts));
 h82 = uicontrol('Style','push','Units','normalized',...
  'Position',[0.4 0.02 0.12 0.07],'Callback', 'uiresume(gcbf);',...
  'String','Done');
 uiwait(hfig8);
 acctol(6) = str2num(get(c84, 'String'));

 close(hfig8);

end

%------------------------------------------------------------------------------
if editsection(8)~=0
% optimization information: opt

 hfig91 = figure('Name','Miser3.2: Optimization information',...
  'Menubar','none','NumberTitle','off');
 message91={' Edit, where necessary, the information below';
  'and then click ''Done''';
  'The values shown are the default values'};
 h91a = uicontrol('Style','text','FontUnits','points','FontSize',12, ...
  'Units','normalized','Position',[0.2 0.82 0.6 0.14],'String',message91);
 String91={'Maximum number of iterations',...
  'Maximum number of function calls',...
  'Gradient accuracy parameter',...
  'Constraint accuracy requirement'};
 String92={'100';'200';'1e-05';'1e-07'};
 for i=1:4
  c91c(i) = uicontrol('Style','text','Units','normalized',...
   'Position',[0.15 0.7-0.06*(i-1) 0.5 0.05],'String',String91{i});
  c91d(i) = uicontrol('Style','edit','Units','normalized',...
   'Position',[0.7 0.7-0.06*(i-1) 0.1 0.05],'String',String92{i});
 end
 h91c=uicontrol('Style','text','Units','normalized',...
  'Position',[0.05 0.4 0.32 0.05],'String','Printout for optimization:');
 c91e=uicontrol('Style','checkbox','Units','normalized',...
  'Position',[0.42 0.43 0.25 0.05],'String','No output');
 c91f=uicontrol('Style','checkbox','Units','normalized',...
  'Position',[0.72 0.43 0.25 0.05],'String','Final output only');
 c91g=uicontrol('Style','checkbox','Units','normalized',...
  'Position',[0.42 0.36 0.25 0.05],'String','Intermediate output','Value',1);
 c91h=uicontrol('Style','checkbox','Units','normalized',...
  'Position',[0.72 0.36 0.25 0.05],'String','Detailed output');
 h91d=uicontrol('Style','text','Units','normalized',...
  'Position',[0.02 0.23 0.3 0.05],'String','Optimization line search:');		
 c91i=uicontrol('Style','checkbox','Units','normalized',...
  'Position',[0.35 0.26 0.4 0.05],'String',...
  'quadratic/cubic polynomial','Value',1);
 c91j=uicontrol('Style','checkbox','Units','normalized',...
  'Position',[0.35 0.2 0.6 0.05],'String',...
  'cubic polynomial (not for numerical derivatives)');		
 h91g = uicontrol('Style','push','Units','normalized',...
  'Position',[0.4 0.05 0.1 0.07],'Callback', 'uiresume(gcbf);','String','Done');
 uiwait(hfig91);

 for i=1:2
  opt(i)= str2num(get(c91d(i), 'String'));
 end
 for i=1:2
  opt(i+3)= str2num(get(c91d(i+2), 'String'));
 end
 opt(3)=0;
 if get(c91f, 'Value')==1
  opt(3)=1;
 elseif get(c91g, 'Value')==1
  opt(3)=2;
 else
  opt(3)=3;
 end
 opt(6)=get(c91j, 'Value');
	
 close(hfig91);

end

%------------------------------------------------------------------------------
if editsection(9)~=0
% numerical derivative information: numder

 hfig93 = figure('Name','Miser3.2: User Derivative information',...
  'Menubar','none','NumberTitle','off');
 message93a={' Edit, if necessary, the information below',...
  ' and then click ''Done'''};
 message93b={' Derivatives of user supplied functions:',...
  '0--analytic (user supplied); 1--numerical'};
 message93c={'f             g0            phi             g             x0'};
 message93d={'x derivative','u derivative','z derivative'};
 h93a = uicontrol('Style','text','FontUnits','points','FontSize',14,...
  'Units','normalized','Position',[0.2 0.8 0.6 0.12],'String',message93a);
 h93b=uicontrol('Style','text','FontUnits','points','FontSize',12,...
  'Units','normalized','Position',[0.2 0.62 0.6 0.1],'String',message93b);
 h93c=uicontrol('Style','text','Units','normalized',...
  'Position',[0.3 0.5 0.45 0.05],'String',message93c);
 message93d={'x derivatives','u derivatives','z derivatives'};
 for i=1:3
  h93d=uicontrol('Style','text','Units','normalized',...
   'Position',[0.1 0.44-(i-1)*0.06 0.16 0.05],'String',message93d(i));
  for j=1:5
   if (i==1 & j~=5) | (i==2 & j==1|j==2|j==4) | i==3
    c93(i,j)=uicontrol('Style','edit','Units','normalized', ...
     'Position',[0.3+(j-1)*0.1 0.44-(i-1)*0.06 0.05 0.05],'String','0');
   end
  end
 end
 h93e = uicontrol('Style','push','Units','normalized',...
  'Position',[0.4 0.2 0.12 0.06],'Callback', 'uiresume(gcbf);','String','Done');
 uiwait(hfig93);

 numder=zeros(3,5);
 for j=1:4
  numder(1,j)= str2num(get(c93(1,j), 'String'));
 end
 for j=1:2
  numder(2,j)= str2num(get(c93(2,j), 'String'));
 end
 numder(2,4)= str2num(get(c93(2,4), 'String'));
 for j=1:5
  numder(3,j)= str2num(get(c93(3,j), 'String'));
 end	
 close(hfig93);

end

%------------------------------------------------------------------------------
if editsection(10)~=0
% file information: kout & filenames

 len=find(fname2=='.')-1;
 fnameshort=fname2(1:len);
 kout=[1 1 0 0];
 hfig94 = figure('Name','Miser3.2: Input, Output and Error Information',...
  'Menubar','none','NumberTitle','off');
 message94=...
  {' Edit, where necessary, the information below and then click ''Done''';
  'The values shown are the default values'};
 h94a = uicontrol('Style','text','FontUnits','points','FontSize',12, ...
  'Units','normalized','Position',[0.05 0.89 0.9 0.1],'String',message94);

 h94b=uicontrol('Style','text','FontUnits','points','FontSize',11,...
  'Units','normalized','Position',[0.05 0.8 0.32 0.07],...
  'String','Error message control:');
 c94a=uicontrol('Style','checkbox','Units','normalized', ...
  'Position',[0.4 0.83 0.2 0.05],'String','Disaster only');
 c94b=uicontrol('Style','checkbox','Units','normalized', ...
  'Position',[0.65 0.83 0.2 0.05],'String','Warnings');
 c94c=uicontrol('Style','checkbox','Units','normalized', ...
  'Position',[0.4 0.78 0.2 0.05],'String','Some debug','Value',1);
 c94d=uicontrol('Style','checkbox','Units','normalized', ...
  'Position',[0.65 0.78 0.2 0.05],'String','All debug');
 h94c=uicontrol('Style','text','Units','normalized',...
  'Position',[0.3 0.69 0.25 0.05],'String','Error file name:');
 c94e=uicontrol('Style','edit','Units','normalized',...
  'Position',[0.55 0.69 0.2 0.05],'String',strcat(fnameshort,'.err'));

 h94d=uicontrol('Style','text','FontUnits','points','FontSize',11, ...
  'Units','normalized','Position',[0.05 0.605 0.2 0.05],...
  'String','Restart file:');
 c94f=uicontrol('Style','checkbox','Units','normalized', ...
  'Position',[0.1 0.55 0.1 0.05],'String','No');
 c94g=uicontrol('Style','checkbox','Units','normalized', ...
  'Position',[0.25 0.55 0.1 0.05],'String','Yes','Value',1);
 h94e=uicontrol('Style','text','Units','normalized',...
  'Position',[0.4 0.55 0.25 0.05],'String','Restart file name');
 c94h=uicontrol('Style','edit','Units','normalized',...
  'Position',[0.65 0.55 0.2 0.05],'String',strcat(fnameshort,'.res'));

 h94f=uicontrol('Style','text','FontUnits','points','FontSize',11, ...
  'Units','normalized','Position',[0.05 0.455 0.15 0.05],...
  'String','Save file:');
 c94i=uicontrol('Style','checkbox','Units','normalized', ...
  'Position',[0.1 0.4 0.1 0.05],'String','No');
 c94j=uicontrol('Style','checkbox','Units','normalized', ...
  'Position',[0.25 0.4 0.1 0.05],'String','Yes','Value',1);
 h94g=uicontrol('Style','text','Units','normalized',...
  'Position',[0.4 0.4 0.25 0.05],'String','Save file name');
 c94k=uicontrol('Style','edit','Units','normalized',...
  'Position',[0.65 0.4 0.2 0.05],'String',strcat(fnameshort,'.sav'));
		
 h94h=uicontrol('Style','text','FontUnits','points','FontSize',11,...
  'Units','normalized','Position',[0.05 0.28 0.28 0.05],...
  'String','Solution file control:');
 c94l=uicontrol('Style','checkbox','Units','normalized', ...
  'Position',[0.35 0.3 0.5 0.05],'String',...
  'solution file does not contain costates','Value',1);
 c94m=uicontrol('Style','checkbox','Units','normalized', ...
  'Position',[0.35 0.24 0.5 0.05],'String',...
  'solution file contains costates');
 h94i=uicontrol('Style','text','Units','normalized',...
  'Position',[0.3 0.16 0.25 0.05],'String','Solution file name:');
 c94n=uicontrol('Style','edit','Units','normalized',...
  'Position',[0.55 0.16 0.2 0.05],'String',strcat(fnameshort,'.sol'));
		
 h94j = uicontrol('Style','push','Units','normalized',...
  'Position',[0.4 0.02 0.12 0.07],'Callback', 'uiresume(gcbf);',...
  'String','Done');
 uiwait(hfig94);
	
 if get(c94a, 'Value')==1
  kout(1)=0;
 elseif get(c94b, 'Value')==1
  kout(1)=1;
 elseif get(c94d, 'Value')==1
  kout(1)=3;
 end
 kout(3)=get(c94m, 'Value');
 f1=get(c94e, 'String');
 f2=get(c94h, 'String');
 kout(2)=get(c94g, 'Value');
 f3=get(c94n,'String');
 kout(4)=get(c94j, 'Value');
 f4=get(c94k, 'String');
 fileinout=char(f1,f2,f3,f4);
	
 close(hfig94);

end

%------------------------------------------------------------------------------
if editsection(11)~=0
% Miscellaneous information: misc
	
 misc=[];
 hfig97 = figure('Name','Miser3.2: Miscellaneous Information',...
  'Menubar','none','NumberTitle','off');
 message97={' Edit, where necessary, the information below';
  'and then click ''Continue''';'The values shown are the default values'}; 
 h97a = uicontrol('Style','text','FontUnits','points','FontSize',12, ...
  'Units','normalized','Position',[0.2 0.84 0.6 0.15],'String',message97);
 h97b=uicontrol('Style','text','Units','normalized',...
  'Position',[0.05 0.76 0.7 0.05],...
  'String','Frequency of user derivative checks (''9999'' means no checks)');
 c97a=uicontrol('Style','edit','Units','normalized',...
  'Position',[0.8 0.76 0.1 0.05],'String','9999');
 h97c=uicontrol('Style','text','Units','normalized',...
  'Position',[0.05 0.7 0.7 0.05],...
  'String','Frequency of monitor and save (''0'' means no saves)');
 c97b=uicontrol('Style','edit','Units','normalized',...
  'Position',[0.8 0.7 0.1 0.05],'String','10');
 h97d=uicontrol('Style','text','Units','normalized',...
  'Position',[0.05 0.64 0.6 0.05],...
  'String','Stiff ode solver for state equations?');
 c94c=uicontrol('Style','checkbox','Units','normalized', ...
  'Position',[0.7 0.64 0.1 0.05],'String','No','Value',1);
 c94d=uicontrol('Style','checkbox','Units','normalized', ...
  'Position',[0.85 0.64 0.1 0.05],'String','Yes');
 h97e=uicontrol('Style','text','Units','normalized',...
  'Position',[0.05 0.58 0.6 0.05],...
  'String','Stiff ode solver for costate equations?');
 c94e=uicontrol('Style','checkbox','Units','normalized', ...
  'Position',[0.7 0.58 0.1 0.05],'String','No','Value',1);
 c94f=uicontrol('Style','checkbox','Units','normalized', ...
  'Position',[0.85 0.58 0.1 0.05],'String','Yes');
 if sysdata(6)	
  h97f=uicontrol('Style','text','Units','normalized',...
   'Position',[0.05 0.5 0.4 0.05],...
   'String','Specify type of absolute value smoothing:');
  c97e=uicontrol('Style','checkbox','Units','normalized', ...
   'Position',[0.15 0.44 0.3 0.05],'String','(g*g+rho^2)/(2*rho)','Value',1);
  c97f=uicontrol('Style','checkbox','Units','normalized', ...
   'Position',[0.5 0.44 0.25 0.05],'String','sqrt(g*g+rho^2)');
 end	
 h97g=uicontrol('Style','text','Units','normalized',...
  'Position',[0.25 0.37 0.4 0.05],...
  'String','Number of parameters set by the user');
 c97g=uicontrol('Style','edit','Units','normalized',...
  'Position',[0.8 0.37 0.1 0.05],'String','0');		
 h97h = uicontrol('Style','push','Units','normalized',...
  'Position',[0.4 0.3 0.12 0.06],'Callback', 'uiresume(gcbf);',...
  'String','Continue');
 uiwait(hfig97);
	
 misc(1)=str2num(get(c97a, 'String'));
 misc(2)=get(c94d, 'Value');
 misc(3)=get(c94f, 'Value');
 misc(4)=str2num(get(c97b, 'String'));
 misc(5)=1;
 if sysdata(6), misc(5)=1+get(c97f, 'Value'); end
 misc(6)=str2num(get(c97g, 'String'));
 npar=misc(6);
	
 if npar > 0
  h97f=uicontrol('Style','text','Units','normalized',...
   'Position',[0.1 0.22 0.7 0.05],...
   'String','Enter the values of the user parameters, 1 per box');
  for j=1:npar
   if j < 9
    h97g(j) = uicontrol('Style','edit','Units','normalized',...
     'Position',[0.05+(j-1)*0.12 0.16 0.1 0.05]);
   else
    h97g(j) = uicontrol('Style','edit','Units','normalized',...
     'Position',[0.05+(j-9)*0.12 0.1 0.1 0.05]);
   end
  end 
  h97e = uicontrol('Style','push','Units','normalized',...
   'Position',[0.4 0.02 0.12 0.07],'Callback', 'uiresume(gcbf);',...
   'String','Done');
  uiwait(hfig97);
		
  for j=1:npar
   misc(6+j)=str2num(get(h97g(j), 'String'));
  end
 end

 close(hfig97);

end

% End of editdata
		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function writefile(fid,title,lcdate,typefi,rsinfo,sysdata,knotdata,...
 condata1,condata2,pardata,consdata1,consdata2,acctol,opt,numder,kout,...
 fileinout,misc)

% Writes the problem input data file

% Write file information
str1='*************************';
fprintf(fid,'%s       File Information      %s\n',str1,str1);
fprintf(fid,' File Title: %s\n',title);
fprintf(fid,' Last Change: %s\n',lcdate);
fprintf(fid,' File Type: %s\n',typefi);
version=' 2.0 ';
fprintf(fid,' MISER3 Version: %s\n',version);
fprintf(fid,' Restart Info: %s\n',rsinfo);

% Write system information
fprintf(fid,'\n%s      System Information     %s\n',str1,str1);
fprintf(fid,' ns =     %4i           (number of states)\n',sysdata(1));
fprintf(fid,' nc =     %4i           (number of controls)\n',sysdata(2));
fprintf(fid,' nz =     %4i           (number of system parameters)\n',...
 sysdata(3));
fprintf(fid,' tstart %16.8e (objective functional integral lower limit)\n',...
 sysdata(4));
fprintf(fid,' tfinal %16.8e (objective functional integral upper limit)\n',...
 sysdata(5));
labch='F';
if sysdata(6)==1,labch='T';end
fprintf(fid,' lab = %c                 (objective or constraints use',labch);
fprintf(fid,' abs smoothing)\n');

% Write knot sets definition
nksets=sysdata(7);
if nksets == 0
 kntype=[];
 knset=[];
else
 kntype=knotdata(:,1);
 knset=knotdata(:,2);
end
fprintf(fid,'\n%s     Knot Sets Definition    %s\n',str1,str1);
fprintf(fid,' nksets = %4i           (number of knot sets)\n',nksets);
fprintf(fid,' knot set #     knot type     number of knots\n');
for i=1:sysdata(7)
 fprintf(fid,' %5i              %1i           %6i\n',i,kntype(i),knset(i));
end
fprintf(fid,'\n%s  Knots for Type 1 Knot Sets %s\n',str1,str1);
if sysdata(7)>0
 knindex=find(kntype==1);
 if length(knindex) > 0
  for i=1:length(knindex)
   fprintf(fid,'knot set #%2i',knindex(i));
   fprintf(fid,'\n%15.7e%15.7e%15.7e%15.7e%15.7e',...
   knotdata(knindex(i),3:2+knset(knindex(i))));
   fprintf(fid,'\n');
  end
 end
end

% Write control definition
nc=sysdata(2);
str2='*****************';
if nc > 0
 kkc=condata1(:,1)';
 koncts=condata1(:,2)';
 kreg=condata1(:,3:11);
 fprintf(fid,'\n%s      Control Definition     %s\n',str1,str1);
 fprintf(fid,' control #   knot set   continuity  regularization type\n');
 for i=1:nc
  fprintf(fid,'    %2i          %2i          %2i',i,kkc(i),koncts(i));
  fprintf(fid,'            %2i %2i %2i  %2i %2i %2i  %2i %2i %2i\n',...
   kreg(i,:));
 end
   
% Write control lower, initial and upper values and regularization weights
 fprintf(fid,'\n*%s Control lower initial upper regularization %s\n',...
  str2,str2);
 for ic=1:nc, kncp(ic)=knset(kkc(ic))+koncts(ic)-1; end
 maxk=max(kncp);
 for i=1:nc
  fprintf(fid,' control #%2i\n',i);
  for j=1:kncp(i)
   fprintf(fid,'%4i  %12.4e%16.8e%12.4e%10.2e%10.2e%10.2e \n',...
    j,condata2(i,j),condata2(i,maxk+j),condata2(i,2*maxk+j),...
    condata2(i,3*maxk+j),condata2(i,4*maxk+j),condata2(i,5*maxk+j));
  end
 end
end

% Write system parameters
nz=sysdata(3);
if nz > 0
 ncp=nc+1;
 fprintf(fid,'\n%s      System Parameters      %s\n',str1,str1);
 fprintf(fid,' No.  lower bound     initial value   upper bound\n');
 for i=1:nz
  fprintf(fid,'%5i     %16.8e%16.8e%16.8e\n',i,pardata(i),pardata(nz+i),...
   pardata(2*nz+i));
 end
end

% Write constraint information
fprintf(fid,'\n%s    Constraint Information   %s\n',str1,str1);
fprintf(fid,'%5i%5i%5i           (Equality constraints of type)\n',...
 consdata1(2:4));
fprintf(fid,'%5i%5i%5i           (Inequality constraints of type)\n',...
 consdata1(5:7));
fprintf(fid,'No.  Type            Characteristic Time  Eps-Tau Smooth');
fprintf(fid,'  Abs Smoothing');
ngu=consdata1(2)+consdata1(5);
ngc=consdata1(3)+consdata1(6);
ng=consdata1(1);
if ng > 0
 contype=consdata2(1:ng);
 taui=ones(1,ng)*sysdata(5);
 index=find(contype==2 | contype==12 | contype==13);
 lindex=length(index);
 index1=find(contype==1 | contype==11);
 i1=length(index1)*(nc+1)+ng;
 i2=i1+lindex;
 i3=i2+lindex;
 tau=consdata2(i1+1:i2);
 taui(index)=tau;
 lget=consdata2(i2+1:i3);
 lgab=consdata2(i3+1:i3+lindex);
 ig=0;iz=ng;
 for j=1:ng
  etm='na ';
  absm='na ';
  if contype(j)==1,descr='lin.eq.u  ';end
  if contype(j)==2
   ig=ig+1;
   descr='eq.canon  ';
   absm='no ';
   if lgab(ig)==1,absm='yes';end   
  end
  if contype(j)==3,descr='eq.system ';end
  if contype(j)==11,descr='lin.ineq.u';end
  if contype(j)==12
   ig=ig+1;
   descr='ineq.canon';
   absm='no ';
   if lgab(ig)==1,absm='yes';end   
  end
  if contype(j)==13
   ig=ig+1;
   descr='ineq.allti';
   absm='no ';
   if lgab(ig)==1,absm='yes';end 
   etm='no ';
   if lget(ig)==1,etm='yes';end
  end
  if contype(j)==14,descr='ineq.sys  ';end
  fprintf(fid,'\n%3i   %2i %s  %16.8e          %s          %s',...
   j,contype(j),descr,taui(j),etm,absm);
  if contype(j)==1|contype(j)==11
   fprintf(fid,'\n%16.8e%16.8e%16.8e%16.8e%16.8e',consdata2(iz+1:iz+nc+1));
   iz=iz+nc+1;
  end
 end
end

% Write accuracy and tolerance
fprintf(fid,'\n\n%s    Accuracy and Tolerance   %s\n',str1,str1);
fprintf(fid,' tolx = %16.8e (state equation tolerance)\n',acctol(1));
fprintf(fid,' tolpsi %16.8e (co-state equation tolerance)\n',acctol(2));
fprintf(fid,' hmax = %16.8e (quadrature interval relative maximum)\n',...
 acctol(3));
fprintf(fid,' reltest%16.8e (relative error test)\n',acctol(4));
fprintf(fid,' epsjts %16.8e (eps of eps-tau algorithm)\n',acctol(5));
fprintf(fid,' taujts %16.8e (tau of eps-tau algorithm)\n',acctol(6));
fprintf(fid,' rhoabs %16.8e (rho of abs smoothing)\n',acctol(7));

% Write optimization selection
fprintf(fid,'\n%s    Optimization Selection   %s\n',str1,str1);
fprintf(fid,' maxite=  %4i           (Maximum Number of Iterations)\n',opt(1));
fprintf(fid,' maxfun=  %4i           (Maximum Number of Function Calls)\n',...
 opt(2));
fprintf(fid,' optprt=  %4i           (Optimization reporting level on',opt(3));
fprintf(fid,' printing)\n');
fprintf(fid,' epsopt=%16.8e (Accuracy for variables and function)\n',opt(4));
fprintf(fid,' epscon=%16.8e (Constraint accuracy)\n',opt(5));
fprintf(fid,' lsearch= %4i           (Line Search - quad/cubic (0);',opt(6));
fprintf(fid,' cubic (1))\n');

% Write user derivatives selection
fprintf(fid,'\n%s    Derivatives of User Supplied Functions   %s\n',str2,str2);
fprintf(fid,' 0--analytic (user supplied) derivatives;  1--numerical');
fprintf(fid,' derivatives\n');
fprintf(fid,'                  f    g0    phi    g    x0\n');
fprintf(fid,' x derivatives    %i     %i     %i     %i      \n',numder(1,1:4));
fprintf(fid,' u derivatives    %i     %i     ',numder(2,1:2));
fprintf(fid,'      %i      \n',numder(2,4));
fprintf(fid,' z derivatives    %i     %i     %i     %i     %i\n',numder(3,1:5));

% Write input, output and error file info
fprintf(fid,'\n%s   Input, Output and Error   %s\n',str1,str1);
fprintf(fid,' error messages:   kout1= %2i  fnout1= %s\n',kout(1),...
 fileinout(1,:));                                    
fprintf(fid,' restart file:     kout2= %2i  fnout2= %s\n',kout(2),...
 fileinout(2,:));                                     
fprintf(fid,' solution file:    kout3= %2i  fnout3= %s\n',kout(3),...
 fileinout(3,:));                                                                          
fprintf(fid,' save file:        kout4= %2i  fnout6= %s\n',kout(4),...
 fileinout(4,:));

% Write miscellaneous selection
fprintf(fid,'\n%s        Miscellaneous        %s\n',str1,str1);
fprintf(fid,' ncheck = %4i           (Frequency of user derivative check)\n',...
 misc(1));
fprintf(fid,' ststiff =%4i           (State ode solver - ode45 (0);',misc(2));
fprintf(fid,' ode15s (stiff) (1))\n');
fprintf(fid,' costiff =%4i           (Costate ode solver - ode45 (0);',misc(3));
fprintf(fid,' ode15s (stiff) (1))\n');
fprintf(fid,' nsave =  %4i           (How often to save)\n',misc(4));
fprintf(fid,' kabs =   %4i           (Which absolute value smoothing)\n',...
 misc(5));
fprintf(fid,' nupar =  %4i           (How many user parameters)\n',misc(6));
if misc(6) > 0
 fprintf(fid,' %15.7e%15.7e%15.7e%15.7e%15.7e\n',misc(7:end));
end

fprintf(fid,'\n%s       End of Data File      %s\n',str1,str1);

% End of writefile
		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [koncts,kpcp,kncp,cspar,clower,cupper,lreg,kreg,wreg0,wreg1,wreg2,h,...
 tk,dtk,tknot,kptk,kpqk,tqd,tau,kpqtau,kpar,Auin,Buin,Aueq,Bueq,ngui,ferr]=...
 setup(ns,nc,nz,tstart,tfinal,nksets,knotdata,condata1,condata2,pardata,...
 consdata1,consdata2,hmax,title,lcdate,typefi,version,rsinfo,fileinout)
%  To set all internal parameters and pointers from user specifications.

% unpack some of the problem variables
if nksets > 0
 knset=knotdata(:,2);
 atk=knotdata(:,3:end);
end
koncts=[];
if nc > 0
 kkc=condata1(:,1)';
 koncts=condata1(:,2)';
end
if nz > 0
 zlwbnd=pardata(:,1);
 zinival=pardata(:,2);
 zupbnd=pardata(:,3);
end
ng=consdata1(1);
ngueq=consdata1(2);
ngceq=consdata1(3);
ngzeq=consdata1(4);
nguin=consdata1(5);
ngcin=consdata1(6);
ngzin=consdata1(7);
ngu=ngueq+nguin;
ngc=ngceq+ngcin;
if ngu > 0
 for i=1:ngu
  alpha(nc*(i-1)+1:i*nc)=consdata2(ng+(nc+1)*(i-1)+1:ng+(nc+1)*(i-1)+nc);
  beta(i)=consdata2(ng+(nc+1)*i);
 end
end
if ngc > 0
 tau=consdata2(ng+ngu*(nc+1)+1:ng+ngu*(nc+1)+ngc);
else
 tau=[];
end

fname1=deblank(fileinout(1,:));
fname2=deblank(fileinout(4,:));
if strcmp(fname1,fname2)==0
 ferr=fopen(fname1,'w');
 fprintf('\n Opening file %s for error reporting',fname1);
 fprintf(ferr,' File Title: %s\n',title);
 d=datenum(now);
 d1=datestr(d,'ddd');d2=datestr(d,'mmm');d3=datestr(d,'dd');
 d4=datestr(d,'HH:MM:SS');d5=datestr(d,'yyyy');
 tdate=char(strcat(d1,{' '},d2,{' '},d3,{' '},d4,{' '},d5));
 fprintf(ferr,'\n Todays Date: %s\n',tdate);
 fprintf(ferr,' File Type:           errors\n');
 fprintf(ferr,' MISER3 Version: %s\n',version);
 fprintf(ferr,' Restart Info: %s\n',rsinfo);
end

%______________________________________________________________________________
% Set up all parameters and variables
bigno=1.0e20;
srepsm=sqrt(eps);

%______________________________________________________________________________
% Set up the array tk which contains all knots for the knot sets 
% in knot set order. 

% set kptkset - pointers to start of eack knot set.
kptkset(1)=1;
for ik=2:nksets
 kptkset(ik)=kptkset(ik-1)+knset(ik-1);
end

%     -- set up tk --
tk=[tstart tfinal];h=[];
for i=1:nksets
 tk(kptkset(i):kptkset(i)+knset(i)-1)=atk(i,1:knset(i));
 h(kptkset(i))=-bigno;
 h(kptkset(i)+1:kptkset(i)+knset(i)-1)=atk(i,2:knset(i))-atk(i,1:knset(i)-1);
end

% set no. of parameters for each control
nct=0;lreg=0;
kpcp=[];kptk=[];kncp=[];kreg=[];wreg0=[];wreg1=[];wreg2=[];
if nc > 0
 for ic=1:nc, kncp(ic)=knset(kkc(ic))+koncts(ic)-1; end
 nct=sum(kncp);
 kptk=kptkset(kkc);
 kpcp(1)=1;
 for i=2:nc
  kpcp(i)=kpcp(i-1)+kncp(i-1);
 end

% set the lower, initial and upper values for the control parameters.
 maxk=max(kncp);
 for i=1:nc
  ncp=kncp(i);
  cspar(kpcp(i):kpcp(i)+ncp-1)=condata2(i,maxk+1:maxk+ncp);
  clower(kpcp(i):kpcp(i)+ncp-1)=condata2(i,1:ncp);
  cupper(kpcp(i):kpcp(i)+ncp-1)=condata2(i,2*maxk+1:2*maxk+ncp);
 end
 nct=length(cspar);

% set the regularization weights for the control parameters.
 kreg=condata1(:,3:11);
 wreg0=zeros(nc,maxk);wreg1=zeros(nc,maxk);wreg2=zeros(nc,maxk);
 for i=1:nc
  if any(kreg(i,:)+1), lreg=1; end;
  ncp=kncp(i);
  wreg0(i,1:ncp)=condata2(i,3*maxk+1:3*maxk+ncp);
  wreg1(i,1:ncp)=condata2(i,4*maxk+1:4*maxk+ncp);
  wreg2(i,1:ncp)=condata2(i,5*maxk+1:5*maxk+ncp);
 end
end

%  set the system parameters lower, initial and upper.
if nz > 0
 cspar(nct+1:nct+nz)=zinival(1:nz);
 clower(nct+1:nct+nz)=zlwbnd(1:nz);
 cupper(nct+1:nct+nz)=zupbnd(1:nz);	
end
nctot=nct+nz;

%______________________________________________________________________________
% set up internal knot set of all knots and characteristic times.
tknot=union(tk,tau);
nknot=length(tknot);
	
%  set up quadrature points.
%  must be a multiple of 4 subintervals between consecutive internal
%  knots. User specifies maximum hmax (relative) for quadratures.
nqp=1;
tqd(1)=tknot(1);
ahmax=hmax*abs(tfinal-tstart);
for ik=1:nknot-1
 dt=tknot(ik+1)-tknot(ik);	
 numstep=fix((abs(dt/ahmax)-0.999)/4)*4+4;
 dtk(ik)=dt/numstep;
 tqd(nqp+1:nqp+numstep)=[tknot(ik)+dtk(ik):dtk(ik):tknot(ik+1)];
 kpqk(ik)=nqp;
 nqp=nqp+numstep;
end	
kpqk(nknot)=nqp;	
	
% set pointers for characteristic times into quadrature points.
% remember that tau(0)=tfinal and kpqtau(0) are not included
kpqtau=[];
for i=1:ngc
 [dum,ia,ib]=intersect(tau(i),tqd);
 kpqtau(i)=ib;
end

% set kpar - internal knot subinterval map to parameters for piecewise
% constant and piecewise linear basis functions.
kpar=[];
for ic=1:nc
 k=0;
 for ik=1:nknot-1
  if(tknot(ik) < tk(kptk(ic)+k) | tknot(ik+1) > tk(kptk(ic)+k+1))
   k=k+1;
  end
  kpar(ic,ik)=kpcp(ic)+k;
 end
end

%______________________________________________________________________________
% set up coefficients of linear constraints of controls only. 
% The inequality constraints of the form A*u+b>=0 have to be transformed
% to inequality constraints of the form Auin*par<=Buin.
%		Auin is nguin*(nknot-1) x nctot if no piecewise linear controls
%		Auin is nguin*nknot x nctot if at least one piecewise linear control
% equality constraints are of the form Aueq*par=Bueq; 
%		Aueq is ngueq*(nknot-1) x nctot if no piecewise linear controls
%		Aueq is ngueq*nknot x nctot if at least one piecewise linear control

Aueq=[];Bueq=[];
Auin=[];Buin=[];
pl=0;
k=0;
if sum(koncts)>0, pl=1;end
if ngu > 0
 numrows=ngu*(nknot-1);
 if pl, numrows=ngu*nknot; end
 Au=zeros(numrows,nctot);
 ja=0;
 for ig=1:ngu
  for ik=1:nknot-1
   k=k+1;
   tt=tknot(ik);
   for ic=1:nc
    if koncts(ic)==0
% piecewise constants
     Au(k,kpar(ic,ik))=alpha(ja+ic);
    else
% piecewise linear continuous.
% find the knots straddling internal knot.
     ind=find(tk(kptk(ic):kptk(ic)+knset(kkc(ic))) > tt);
     j=ind(1)-1;
     ticj=tk(j);
     ticjp=tk(j+1);
     Au(k,kpar(ic,ik))=alpha(ja+ic)*(ticjp-tt)/(ticjp-ticj);
     Au(k,kpar(ic,ik)+1)=alpha(ja+ic)*(tt-ticj)/(ticjp-ticj);
    end
   end
   Bu(k)=beta(ig);
  end		
  if pl
% add last knot (tfinal) as at least one piecewise linear control.
   k=k+1;
   Bu(k)=beta(ig);
   for ic=1:nc
    Au(k,kncp(ic)+kpcp(ic)-1)=alpha(ja+ic);
   end
  end
  ja=ja+nc;		
 end
 ngui=k; 
 k=(ngui/ngu);
else
 ngui=0;
end
% Form the matrices for equality and inequality constraints
ngueqi=ngueq*k;
nguini=nguin*k;
if ngueq > 0
 Aueq=Au(1:ngueqi,:);
 Bueq=Bu(1:ngueqi)';
end
if nguin > 0
 Auin=-Au(ngueqi+1:ngui,:);
 Buin=Bu(ngueqi+1:ngui)';
end

% End of setup

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function kountsav=ocsave(ktype,iter,epsjt,taujt,rhoab,lab,let,kpcp,kncp,...
 cspar,title,lcdate,typefi,version,rsinfo,sysdata,knotdata,condata1,...
 condata2,pardata,consdata1,consdata2,acctol,opt,numder,kout,fileinout,misc)
% This routine saves the current parameters to the save file, by 
% appending two heading lines, a data input file and two tail lines.

global Once 

% ktype=0 save when iter is divisible by nsave.
%      =1 save at end of a sequential optimization.
%      =2 save at end of run, same as restart file.
%______________________________________________________________________________
nsave=misc(4);

if kout(4)*nsave > 0 & (ktype > 0 | rem(iter,nsave)==0)

 rw='a';
 if Once==0, rw='w'; end
 fname=deblank(fileinout(4,:));
 fidsav=fopen(fname,rw);
 if Once==0
  fprintf('\n Opening save file %s for writing',fname);
  fprintf(fidsav,' File Title: %s\n',title);
  fprintf(fidsav,' Last Change: %s\n',lcdate);
  fprintf(fidsav,' File Type:       Save file\n');
  fprintf(fidsav,' MISER3 Version: %s\n',version);
  fprintf(fidsav,' Restart Info: %s\n',rsinfo);
 end
%______________________________________________________________________________
% Put current parameters into the initial arrays.
 maxk=max(kncp);
 nct=sum(kncp);
 nz=sysdata(3);
 for i=1:sysdata(2)
  ncp=kncp(i);
  condata2(i,maxk+1:maxk+ncp)=cspar(kpcp(i):kpcp(i)+ncp-1);
 end
 if nz > 0
  zinival=cspar(nct+1:nct+nz);
  pardata(:,2)=zinival';
 end
%______________________________________________________________________________
 type='iteration ';
 if ktype==1 & lab, type='end rhoabs'; end
 if ktype==1 & let, type='end epstau'; end
 if ktype==2, type='final'; end
 fprintf(fidsav,'\n BEGIN SAVE SET, count=%4i, TYPE of SAVE %2i  %s',...
  Once,ktype,type);
 kountsav=0;
 if ktype==0
  fprintf(fidsav,'\n Iteration number is %4i',iter);
 elseif ktype==1
  kountsav=Once;
  fprintf(fidsav,'\n ih,epsjt,taujt,rhoab are %3i%16.8e%16.8e%16.8e%16.8e',...
   iter,epsjt,taujt,rhoab);
  fprintf(fidsav,'\n');
 else
  fprintf(fidsav,'\n FINAL SAVE SET - SAME AS *.RES \n');
 end
		
 writefile(fidsav,title,lcdate,typefi,rsinfo,sysdata,knotdata,condata1,...
  condata2,pardata,consdata1,consdata2,acctol,opt,numder,kout,...
  fileinout,misc);
		
 fprintf(fidsav,' END OF SAVE SET');
 str1='$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$';
 fprintf(fidsav,'\n %s%s\n',str1,str1);
 fclose(fidsav);	
 Once=Once+1;
end

% End of ocsave

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function choice=failurechoice(nfail,let,lab)
% Function failurechoice determines the action to be taken on failure of FMINCON
% if choice = 0, stop
%           = 1, no sequential opt-restart the opt
%           = 2, continue with sequential opt

hfig = figure('Name','Miser3.2: Optimal Control Software',...
 'Menubar', 'none','NumberTitle','off'); 

str1=sprintf('parameter exitflag equal to %3i',nfail);
message1={'Optimization routine FMINCON has failed with failure';str1};
h1 = uicontrol('Style','text','FontUnits','points','FontSize',12,...
 'Units','normalized','Position',[0.1 0.7 0.8 0.1],'String',message1);
str2={'Stop the optimisation','Restart the optimisation',...
 'Continue to the next sequential optimization'};
		
if (let | lab)
 message2={'You are executing a sequential optimization.';
  'Please make one of the choices below then click ''Done''.'};
 h2=uicontrol('Style','text','Units','normalized',...
  'Position',[0.1 0.55 0.8 0.09],'String',message2);
 for i=1:3
  c(i) = uicontrol('Style','checkbox','Units','normalized',...
   'Position',[0.3 0.48-0.07*(i-1) 0.5 0.05],'String',str2{i});
 end
else
 message3={'Please make one of the choices below then click ''Done'''};
 h3=uicontrol('Style','text','Units','normalized',...
  'Position',[0.1 0.48 0.8 0.05],'String',message3);
 for i=1:2
  c(i) = uicontrol('Style','checkbox','Units','normalized',...
   'Position',[0.3 0.41-0.07*(i-1) 0.5 0.05],'String',str2{i});
 end
end

h4 = uicontrol('Style','push','Units','normalized',...
	'Position',[0.44 0.2 0.12 0.07],...
	'Callback', 'uiresume(gcbf);','String','Done');

uiwait(hfig)

cnum(1)=get(c(1), 'Value');
cnum(2)=get(c(2), 'Value');
cnum(3)=0;
if (let | lab), cnum(3)=get(c(3), 'Value');end
if cnum(3)==0
 choice=cnum(2);
else
 choice=2;
end

close(hfig)
		
% End of failurechoice
		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function kresta=restartchoice
% Function restratchoice determines whether to double or triple the control
% parameters in the restart file
% if kresta <= 1, use the same control parameters
%            = 2, double the control parameters
%            = 3, triple the control parameters

hfig = figure('Name','Miser3.2: Optimal Control Software',...
 'Menubar', 'none','NumberTitle','off'); 

h1 = uicontrol('Style','text','FontUnits','points','FontSize',14,...
 'Units','normalized','Position',[0.1 0.85 0.8 0.06],'String',...
 'MISER3 HAS FINISHED');
message1={'Choose the number of control parameters for the restart file';
 '(Default is same number of control parameters as before)'};
h2 = uicontrol('Style','text','FontUnits','points','FontSize',12,...
 'Units','normalized','Position',[0.1 0.68 0.8 0.11],'String',message1);
str2={'Same number of control parameters as before',...
 'Double the number of control parameters',...
 'Triple the number of control parameters'};
c(1) = uicontrol('Style','checkbox','Units','normalized',...
 'Position',[0.2 0.58 0.7 0.05],'String',str2{1},'Value',1);
c(2) = uicontrol('Style','checkbox','Units','normalized',...
 'Position',[0.2 0.5 0.7 0.05],'String',str2{2});
c(3) = uicontrol('Style','checkbox','Units','normalized',...
 'Position',[0.2 0.42 0.7 0.05],'String',str2{3});

h4 = uicontrol('Style','push','Units','normalized',...
 'Position',[0.44 0.28 0.12 0.07],'Callback', 'uiresume(gcbf);',...
 'String','Done');

uiwait(hfig)

kresta=1;
cnum(2)=get(c(2), 'Value');
cnum(3)=get(c(3), 'Value');
if cnum(2)==1
 kresta=2;
elseif cnum(3)==1
 kresta=3;
end

close(hfig)
		
% End of restartchoice
		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ocprint(cspar,f,f0,g,gm,ih,histf,histg,histabs,histeps,histtau,...
 histpar,nfailhist,iterhist,nfail,iterno,funcno,kount,Aueq,Bueq,...
 Auin,Buin,kncp,tknot,lget,kpqk,tqd,title,lcdate,version,...
 rsinfo,sysdata,condata1,consdata1,opt,kout,fileinout,misc,lreg)
% This is the subroutine to output the results for MISER3.

global xqd uqd psiqd

nct=sum(kncp);
nctot=nct+sysdata(3);
if nct > 0, koncts=condata1(:,2); end
nknot=length(tknot);

fname=deblank(fileinout(3,:));
fidsol=fopen(fname,'w');
fprintf('\n  Writing solution file:  %s \n',fname);
fprintf(fidsol,' File Title: %s\n',title);
d=datenum(now);
d1=datestr(d,'ddd');d2=datestr(d,'mmm');d3=datestr(d,'dd');
d4=datestr(d,'HH:MM:SS');d5=datestr(d,'yyyy');
tdate=char(strcat(d1,{' '},d2,{' '},d3,{' '},d4,{' '},d5));
fprintf(fidsol,'\n Todays Date: %s\n',tdate);
fprintf(fidsol,' File Type:        Solution file\n');
fprintf(fidsol,' MISER3 Version: %s\n',version);
fprintf(fidsol,' Restart Info: %s\n',rsinfo);
if misc(6) > 0
 fprintf(fidsol,'\n User parameters for this run are:');
 upar=misc(7:end);
 npar=misc(6);
 nnp=floor((npar-1)/5)+1;
 jmin=1;
 for i=1:nnp
  jmax=min(jmin+4,npar);
  fprintf(fidsol,'\n %15.7e%15.7e%15.7e%15.7e%15.7e',upar(jmin:jmax));
  jmin=jmin+5;
 end
 fprintf(fidsol,'\n');
end

%______________________________________________________________________________
if ih==1
 if (iterno < opt(1) & funcno < opt(2) & nfail > 0)
  fprintf(fidsol,'\n The optimization has finished with no problems in');
  fprintf(fidsol,'%4i iterations \n',iterno);
 end
 if nfail <= 0
  fprintf(fidsol,'\n The optimization routine has failed with exitflag =');
  fprintf(fidsol,'%4i in,%4i iterations \n',nfail,iterno);
 end
 if iterno >= opt(1)
  fprintf(fidsol,'\n The maximum number of iterations was computed');
  fprintf(fidsol,' (%4i).',opt(1));
  fprintf(fidsol,'\n This is the solution at the end of these iterations.\n');
 end	
 if funcno >= opt(2)
  fprintf(fidsol,'\n The maximum number of function calls was computed');
  fprintf(fidsol,' (%4i).',opt(2));
  fprintf(fidsol,'\n This is the solution at the end of these function calls.\n');
 end	
else
 fprintf(fidsol,'\n This is a history of the sequential optimization progress.');
 fprintf(fidsol,'\n no - count number of saveset in save file,');
 fprintf(fidsol,'\n nf - nfail,  it - iterations, ');
 fprintf(fidsol,'\n abs - parameter for absolute value smoothing');
 fprintf(fidsol,' sequential optimization,');
 fprintf(fidsol,'\n eps tau - parameters for smoothing max(o,g)');
 fprintf(fidsol,' sequential optimization,');
 fprintf(fidsol,'\n d param - maximum parameter change (over all');
 fprintf(fidsol,' parameters) since last iteration.\n');
 fprintf(fidsol,'\n  no  nf    it    abs      eps      tau     objective');
 fprintf(fidsol,'        d param');
 itsum=sum(iterhist);
 for i=1:ih
  hdp=max(abs(histpar(i+1,:)-histpar(i,:)));
  fprintf(fidsol,'\n %3i %3i %5i %8.1e %8.1e %8.1e %14.6e %14.6e',...
   kount(i),nfailhist(i),iterhist(i),histabs(i),histeps(i),...
   histtau(i),histf(i),hdp);
 end
 fprintf(fidsol,'\n Total iterations: >=%5i\n',itsum);
	
 if  consdata1(1) > 0
  fprintf(fidsol,'\n History of constraint values at end of each');
  fprintf(fidsol,' sequential optimization step.');
  for i=1:ih
   fprintf(fidsol,'\n %10.2e%10.2e%10.2e%10.2e%10.2e%10.2e%10.2e%10.2e',...
    histg(i,:));
  end
 end		
 fprintf(fidsol,'\n\n Remember the parameters at the end of every');
 fprintf(fidsol,' sequential optimization are\n saved in the save');
 fprintf(fidsol,' file - if you have asked for it.\n Here follows');
 fprintf(fidsol,' the final solution:\n');
end

fprintf(fidsol,'\n The minimum of the OBJECTIVE function is: %16.8e \n',f);
if lreg
 fprintf(fidsol,'\n Value of objective without regularisation: %16.8e \n',f0);
end

ngu=consdata1(2)+consdata1(5);
if ngu > 0
 convaleq=[];convalin=[];
 fprintf(fidsol,'\n Values of linear control only constraints at internal'); 
 fprintf(fidsol,' knots follows. Equality constraints first:');
 if consdata1(2) > 0
  convaleq=Aueq(:,1:nct)*cspar(1:nct)'+Bueq;
 end	
 if consdata1(5) > 0
  convalin=-Auin(:,1:nct)*cspar(1:nct)'+Buin;
 end	
 conval=[convaleq; convalin];
 ip=length(conval)/ngu;
 constrval=reshape(conval,ip,ngu);
 nng=floor((ngu-1)/5)+1;
 jmin=0;
 for j=1:nng
  jmax=min(jmin+4,ngu-1);
  fprintf(fidsol,'\n   Time  Constraint#%2i Constraint#%2i Constraint#%2i Constraint#%2i Constraint#%2i',...
   [jmin+1:jmax+1]);
  for i=1:ip
   fprintf(fidsol,'\n%8.4f%14.6e%14.6e%14.6e%14.6e%14.6e',...
    tknot(i),constrval(i,jmin+1:jmax+1));
  end
  jmin=jmin+5;
 end
end

ngc=consdata1(3)+consdata1(6);
if ngc > 0
 fprintf(fidsol,'\n\n The values of canonical constraints are: ');
 for ig=1:ngc
  if lget(ig)
   fprintf(fidsol,'\n%3i Value of Geps is %12.4e, value of G is %12.4e',...
    ig,g(ig),gm(ig));
  else
    fprintf(fidsol,'\n%3i Value of G is %12.4e',ig,g(ig));
  end
 end
end

ngz=consdata1(4)+consdata1(7);
ng=consdata1(1);	
if ngz > 0
 fprintf(fidsol,'\n The values of system parameter constraints are:');
 for ig=1:ngz
  fprintf(fidsol,'\n%3i Value of G is %12.4e',ig,g(ng-ngc+ig));
 end
end

nc=sysdata(2);
if nc > 0
 fprintf(fidsol,'\n\n The OPTIMAL CONTROLS are given at points called the');
 fprintf(fidsol,' the internal knots. The set\n of internal knots is');
 fprintf(fidsol,' the union of all knot sets and characteristic times.');
 nnc=floor((nc-1)/5)+1;
 jmin=1;
 cts={'  pw constant','  pw linear'};
 for i=1:nnc
  fprintf(fidsol,'\n\n');
  jmax=min(jmin+4,nc);
  fprintf(fidsol,'\n   Time    Control #%2i   Control #%2i   Control #%2i   Control #%2i   Control #%2i',...
   [jmin:jmax]);
  fprintf(fidsol,'\n        ');
  for j=jmin:jmax
   fprintf(fidsol,' %s',char(cts(koncts(j)+1)));
  end
  for j=1:nknot
   fprintf(fidsol,'\n %8.4f%12.5f  %12.5f  %12.5f  %12.5f  %12.5f',...
    tqd(kpqk(j)),uqd(jmin:jmax,kpqk(j)));
  end
  jmin=jmin+5;
 end
end

nz=sysdata(3);
if nz > 0
 fprintf(fidsol,'\n\n OPTIMAL SYSTEM PARAMETERS are:\n');
 fprintf(fidsol,' %13.5e',cspar(nct+1:nct+nz));
end
	
%_____________________________________________________________________________
ns=sysdata(1);
fprintf(fidsol,'\n\n The STATE VARIABLES vs the internal knot set are: ');
nns=floor((ns-1)/5)+1;
jmin=1;
for i=1:nns
 jmax=min(jmin+4,ns);
 fprintf(fidsol,'\n   Time     State #%2i     State #%2i     State #%2i     State #%2i     State #%2i',...
  [jmin:jmax]);
 for j=1:nknot
  fprintf(fidsol,'\n %8.4f%12.5f  %12.5f  %12.5f  %12.5f  %12.5f',...
   tqd(kpqk(j)),xqd(jmin:jmax,kpqk(j)));
 end
 fprintf(fidsol,'\n');
 jmin=jmin+5;
end
	
if kout(3) > 0
 fprintf(fidsol,'\n\n\n\n The Costate variables vs the internal knot set are: ');
 jmin=1;
 for ig=0:ngc
  if ig==0
   fprintf(fidsol,'\n\n Costate variables for Objective function:');
  else
   fprintf(fidsol,'\n\n Costate variables for Constraint number%3i',ig);
  end
  nns=floor((ns-1)/5)+1;
  for is=1:nns
   jmax=min(jmin+4,(ig+1)*ns);
   fprintf(fidsol,'\n   Time   Costate #%2i   Costate #%2i   Costate #%2i   Costate #%2i   Costate #%2i',...
    [jmin:jmax]);
   for j=1:nknot
    fprintf(fidsol,'\n %8.4f%12.5f  %12.5f  %12.5f  %12.5f  %12.5f',...
     tqd(kpqk(j)),psiqd(jmin:jmax,kpqk(j)));
   end
   jmin=jmax+1;
  end
 end
end

fclose(fidsol);
fprintf('\n Finished writing solution file. \n');	

% End of ocprint

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ocrestart(kresta,kpcp,kncp,cspar,title,lcdate,typefi,rsinfo,...
 sysdata,knotdata,condata1,condata2,pardata,consdata1,consdata2,...
 acctol,opt,numder,kout,fileinout,misc)

% This routine saves the current parameters for restart.
% The option to halve or trisect the knot spacing is also included.

nct=sum(kncp);
nc=sysdata(2);
nz=sysdata(3);
nksets=sysdata(7);
if nksets > 0
 knset=knotdata(:,2);
 atk=knotdata(:,3:end);
end

if nc > 0
 koncts=condata1(:,2);
 kreg=condata1(:,3:11);
 maxk=max(kncp);
 lwbnd=zeros(nc,maxk);
 inival=zeros(nc,maxk);
 upbnd=zeros(nc,maxk);
 wreg0=zeros(nc,maxk);wreg1=zeros(nc,maxk);wreg2=zeros(nc,maxk);
 for ic=1:nc
  ncp=kncp(ic);
  lwbnd(ic,1:ncp)=condata2(ic,1:ncp);
  inival(ic,1:ncp)=cspar(kpcp(ic):kpcp(ic)+ncp-1);
  upbnd(ic,1:ncp)=condata2(ic,2*maxk+1:2*maxk+ncp);
  wreg0(ic,1:ncp)=condata2(ic,3*maxk+1:3*maxk+ncp);
  wreg1(ic,1:ncp)=condata2(ic,4*maxk+1:4*maxk+ncp);
  wreg2(ic,1:ncp)=condata2(ic,5*maxk+1:5*maxk+ncp);
 end
end

if nz > 0
 zinival=cspar(nct+1:nct+nz);
 pardata(:,2)=zinival';
end

% Do nothing if no doubling of knot sets.

if kresta==2
% Doubling of parameters on controls by doubling knot sets.
 for ik=1:nksets
  for k=knset(ik):-1:2
   atk(ik,k+k-1)=atk(ik,k);
   atk(ik,k+k-2)=0.5*(atk(ik,k)+atk(ik,k-1));
  end		
  knset(ik)=2*knset(ik)-1;
 end	
 for ic=1:nc
  if koncts(ic)==0
   for ip=kncp(ic):-1:1
    iz=ip+ip;
    lwbnd(ic,iz)=lwbnd(ic,ip);
    lwbnd(ic,iz-1)=lwbnd(ic,ip);
    upbnd(ic,iz)=upbnd(ic,ip);
    upbnd(ic,iz-1)=upbnd(ic,ip);
    inival(ic,iz)=cspar(kpcp(ic)-1+ip);
    inival(ic,iz-1)=inival(ic,iz);
    wreg0(ic,iz)=wreg0(ic,ip);
    wreg0(ic,iz-1)=wreg0(ic,ip);
    wreg1(ic,iz)=wreg1(ic,ip);
    wreg1(ic,iz-1)=wreg1(ic,ip);
    wreg2(ic,iz)=wreg2(ic,ip);
    wreg2(ic,iz-1)=wreg2(ic,ip);
   end
   kncp(ic)=2*kncp(ic);
  else
   for ip=kncp(ic):-1:2
    iz=ip+ip-1;
    lwbnd(ic,iz)=lwbnd(ic,ip);
    lwbnd(ic,iz-1)=0.5*(lwbnd(ic,ip)+lwbnd(ic,ip-1));
    upbnd(ic,iz)=upbnd(ic,ip);
    upbnd(ic,iz-1)=0.5*(upbnd(ic,ip)+upbnd(ic,ip-1));
    inival(ic,iz)=cspar(kpcp(ic)-1+ip);
    inival(ic,iz-1)=0.5*(inival(ic,iz)+cspar(kpcp(ic)-2+ip));
    wreg0(ic,iz)=wreg0(ic,ip);
    wreg0(ic,iz-1)=0.5*(wreg0(ic,ip)+wreg0(ic,ip-1));
    wreg1(ic,iz)=wreg1(ic,ip);
    wreg1(ic,iz-1)=0.5*(wreg1(ic,ip)+wreg1(ic,ip-1));
    wreg2(ic,iz)=wreg2(ic,ip);
    wreg2(ic,iz-1)=0.5*(wreg2(ic,ip)+wreg2(ic,ip-1));
   end			
   inival(ic,1)=cspar(kpcp(ic));
   kncp(ic)=2*kncp(ic)-1;
  end
 end	
end

if kresta==3
% Tripling of parameters on controls by tripling knot sets.
 for ik=1:nksets
  for k=knset(ik):-1:2
   iz=k+k+k-2;
   atk(ik,iz)=atk(ik,k);
   atk(ik,iz-1)=(2*atk(ik,k)+atk(ik,k-1))/3;
   atk(ik,iz-2)=(atk(ik,k)+2*atk(ik,k-1))/3;
  end
  knset(ik)=3*knset(ik)-2;
 end	
 for ic=1:nc
  if koncts(ic)==0
   for ip=kncp(ic):-1:1
    iz=ip+ip+ip;
    lwbnd(ic,iz)=lwbnd(ic,ip);
    lwbnd(ic,iz-1)=lwbnd(ic,ip);
    lwbnd(ic,iz-2)=lwbnd(ic,ip);
    upbnd(ic,iz)=upbnd(ic,ip);
    upbnd(ic,iz-1)=upbnd(ic,ip);
    upbnd(ic,iz-2)=upbnd(ic,ip);
    inival(ic,iz)=cspar(kpcp(ic)-1+ip);
    inival(ic,iz-1)=inival(ic,iz);
    inival(ic,iz-2)=inival(ic,iz);
    wreg0(ic,iz)=wreg0(ic,ip);
    wreg0(ic,iz-1)=wreg0(ic,ip);
    wreg0(ic,iz-2)=wreg0(ic,ip);
    wreg1(ic,iz)=wreg1(ic,ip);
    wreg1(ic,iz-1)=wreg1(ic,ip);
    wreg1(ic,iz-2)=wreg1(ic,ip);
    wreg2(ic,iz)=wreg2(ic,ip);
    wreg2(ic,iz-1)=wreg2(ic,ip);
    wreg2(ic,iz-2)=wreg2(ic,ip);
   end
   kncp(ic)=3*kncp(ic);
  else
   for ip=kncp(ic):-1:2
    iz=ip+ip+ip-2;
    lwbnd(ic,iz)=lwbnd(ic,ip);
    lwbnd(ic,iz-1)=(2*lwbnd(ic,ip)+lwbnd(ic,ip-1))/3;
    lwbnd(ic,iz-2)=(lwbnd(ic,ip)+2*lwbnd(ic,ip-1))/3;
    upbnd(ic,iz)=upbnd(ic,ip);
    upbnd(ic,iz-1)=(2*upbnd(ic,ip)+upbnd(ic,ip-1))/3;
    upbnd(ic,iz-2)=(upbnd(ic,ip)+2*upbnd(ic,ip-1))/3;
    inival(ic,iz)=cspar(kpcp(ic)-1+ip);
    inival(ic,iz-1)=(2*inival(ic,iz)+cspar(kpcp(ic)-2+ip))/3;
    inival(ic,iz-2)=(inival(ic,iz)+2*cspar(kpcp(ic)-2+ip))/3;
    wreg0(ic,iz)=wreg0(ic,ip);
    wreg0(ic,iz-1)=(2*wreg0(ic,ip)+wreg0(ic,ip-1))/3;
    wreg0(ic,iz-2)=(wreg0(ic,ip)+2*wreg0(ic,ip-1))/3;
    wreg1(ic,iz)=wreg1(ic,ip);
    wreg1(ic,iz-1)=(2*wreg1(ic,ip)+wreg1(ic,ip-1))/3;
    wreg1(ic,iz-2)=(wreg1(ic,ip)+2*wreg1(ic,ip-1))/3;
    wreg2(ic,iz)=wreg2(ic,ip);
    wreg2(ic,iz-1)=(2*wreg2(ic,ip)+wreg2(ic,ip-1))/3;
    wreg2(ic,iz-2)=(wreg2(ic,ip)+2*wreg2(ic,ip-1))/3;
   end
   inival(ic,1)=cspar(kpcp(ic));
   kncp(ic)=3*kncp(ic)-2;
  end
 end
end
		
% decrease tolerances.
if kresta > 1
 acctol(1)=acctol(1)*0.1;
 acctol(2)=acctol(2)*0.1;
 opt(4)=opt(4)*0.1;
 opt(5)=opt(5)*0.1;
end

% Put the knot data and control values into the arrays for writing
if nksets > 0
 knotdata(:,2)=knset;
 knotdata=[knotdata(:,1) knset zeros(nksets,max(knset))];
 for i=1:nksets
  knotdata(i,3:knset(i)+2)=atk(i,1:knset(i));
 end
end

if nc > 0
 maxk=max(kncp);
 condata2=zeros(nc,6*maxk);
 for i=1:nc
  ncp=kncp(i);
  condata2(i,1:ncp)=lwbnd(i,1:ncp);
  condata2(i,maxk+1:maxk+ncp)=inival(i,1:ncp);
  condata2(i,2*maxk+1:2*maxk+ncp)=upbnd(i,1:ncp);
  condata2(i,3*maxk+1:3*maxk+ncp)=wreg0(i,1:ncp);
  condata2(i,4*maxk+1:4*maxk+ncp)=wreg1(i,1:ncp);
  condata2(i,5*maxk+1:5*maxk+ncp)=wreg2(i,1:ncp);
 end
end

fname=deblank(fileinout(2,:));
fidre=fopen(fname,'w');
fprintf('\n Writing restart file %s ...\n',fname);
jedit=sscanf(rsinfo(5:8),'%i');
jrest=sscanf(rsinfo(18:21),'%i');
jrest=jrest+1;
rsinfo=sprintf('EDIT%4i, RESTART%4i                                       ',...
 jedit,jrest);	
writefile(fidre,title,lcdate,typefi,rsinfo,sysdata,knotdata,condata1,...
 condata2,pardata,consdata1,consdata2,acctol,opt,numder,kout,fileinout,misc);
fclose(fidre);
if kresta > 1
 kncp=floor(kncp/kresta);
end
jrest=jrest-1;
rsinfo=sprintf('EDIT%4i, RESTART%4i                                       ',...
 jedit,jrest);
fprintf('\n Finished writing restart file.\n');
fprintf('\n');

% End of ocrestart

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotsoln(nc,tqd,tknot,kpqk,koncts)
% plotsoln produces plots of the states and controls

global xqd uqd

% Plot the states against time
if nc > 0, subplot(2,1,1); end
plot(tqd,xqd);
hold on;
[ns dum]=size(xqd);
for i=1:ns, plot(tknot,xqd(i,kpqk),'ko'); end
xlabel('Time (scaled to [0,1] if final time is free)');
ylabel('States');
hold off;

% Plot the controls against time
if nc > 0
 nknot=length(tknot);
 nqp=length(tqd);
 isknot=ismember(tqd,tknot(2:nknot-1));
 i=0;
 for iq=1:nqp
  if isknot(iq)==1
   i=i+1;
   tplot(i)=tqd(iq);
   for ic=1:nc
    if koncts(ic)==0
     uplot(ic,i)=uqd(ic,iq-1);
    else
     uplot(ic,i)=uqd(ic,iq);
    end
   end			
  end
  i=i+1;
  tplot(i)=tqd(iq);
  uplot(:,i)=uqd(:,iq);
 end
 subplot(2,1,2);
 plot(tplot,uplot);
 hold on;
 for i=1:nc, plot(tknot,uqd(i,kpqk),'ko'); end
 xlabel('Time (scaled to [0,1] if final time is free)');
 ylabel('Controls');
 hold off; 
 
 figure(2)
 plot(xqd(1,:),xqd(2,:),'-r');
 hold on;
 plot(500,500,'^r',0,0,'or');
theta=0:0.1:2*pi+0.1;
x0=[452 88 358 497 524;%x
    436 120 78 306 155;%y
    34 26 34 45 53];%r
pho=x0(3,:)-8;

lx=[162,198,380,244,162];
ly=[80,60,185,205,80];
[rlx,rly,k,rb]=rf(lx,ly);

lx1=[127,210,351,234,127];
ly1=[274,245,356,391,274];
[rlx1,rly1,k1,rb1]=rf(lx1,ly1);

hold on
for i=1:size(x0,2)
    x(i,:)=x0(1,i)+x0(3,i)*cos(theta);
    y(i,:)=x0(2,i)+x0(3,i)*sin(theta);
    plot(x(i,:),y(i,:),'k--','HandleVisibility','off')
    rx(i,:)=x0(1,i)+pho(i)*cos(theta);
    ry(i,:)=x0(2,i)+pho(i)*sin(theta);
    patch(transpose(rx(i,:)),transpose(ry(i,:)),'b','HandleVisibility','off')
    %fill(x(i,:),y(i,:));
end
plot(rlx,rly,'k--','HandleVisibility','off')
patch(lx,ly,'b','HandleVisibility','off')
plot(rlx1,rly1,'k--','HandleVisibility','off')
patch(lx1,ly1,'b','HandleVisibility','off')

save('x.mat','xqd');
save('u.mat','uqd');
end

% End of plotsoln

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
