function [xsol,t,J]=v10solver(ic,par,tf,l,varargin)
if nargin==5
    tol=varargin{1};
else
    tol=1e-6;
end
atol=1e-8*ones(1,24);
atol([2,4,6])=1e-16;

if l==1
    options = odeset('RelTol',tol,'AbsTol',atol,'NonNegative', 1:length(ic),...
        'BDF','on','MaxOrder',2,'Stats','off','InitialStep',0.5);
    [t,xout] = ode15s(@prv10,[0 tf],ic,options,par);
    
    xtemp =xout(end,:);
    funb=@(x) prv10(0,x,par);
    opt=optimset('Display','off','MaxFunEvals',1e4,'MaxIter',1e4,'TolFun',...
        1e-6,'TolX',1e-6);
    %[xsol,FVAL,EXITFLAG,OUTPUT,J]=fsolve(funb,xtemp,opt);
    return
elseif l==2
    options = odeset('RelTol',tol,'AbsTol',atol,'NonNegative', 1:length(ic),...
        'BDF','on','MaxOrder',2,'Stats','off','InitialStep',0.5,'Jacobian',@v10sjac); %,'OutputSel',[21,22],'OutputFcn',@odeplot
    [t,xout] = ode15s(@prv10,[0 tf],ic,options,par);
    xtemp =xout(end,:);
    %J=0;
    funb=@(x) fsolprv10(x,par);
    opt=optimset('Display','off','MaxFunEvals',1e4,'MaxIter',1e4,'TolFun',...
        1e-6,'TolX',1e-6,'Jacobian','on'); %
    [xsol,FVAL,EXITFLAG,OUTPUT,J]=fsolve(funb,xtemp,opt);
    return
elseif l==3
    funb=@(x) prv10(0,x,par);
    opt=optimset('Display','iter','MaxFunEvals',1e4,'MaxIter',1e4,'TolFun',...
        1e-6,'TolX',1e-6);
    [xsol,~,~,~,J]=fsolve(funb,ic,opt);
    t=J;
    return
elseif l==4
    options = odeset('RelTol',tol,'AbsTol',atol,'NonNegative', 1:length(ic),...
        'BDF','on','MaxOrder',2,'Stats','off','InitialStep',0.5,'Jacobian',@v10sjac); %,'OutputSel',[21,22],'OutputFcn',@odeplot
    [t,xout] = ode15s(@prv10,[0 tf],ic,options,par);
    xsol=xout(end,:);
    return
elseif l==5
    data.p =par;
    t0 = 0.0;
    y0 =ic(:);
    %tend=time;
    atol=1e-14*ones(24,1);
    atol([2,3,12,9,13,14,16,17,18,19,22,23,24])=1e-16*ones(13,1);
    options = CVodeSetOptions('UserData', data,...
        'RelTol',1e-4,...
        'AbsTol',atol,'MaxNumSteps',1e8,'MaxStep',1e7,...
        'LinearSolver','Dense','JacobianFn',@jv10,'InitialStep',0.5); %
    
    CVodeInit(@rhs, 'BDF', 'Newton', t0, y0, options);
    tout=1:100:tf;
    [status,t,y] = CVode(tout,'Normal');
    CVodeFree;
    y=y(:,end);
    checkneg1=find(y<0, 1);
    if ~isempty(checkneg1)
        warning('negativesolution')
    end
    funb=@(x) fsolprv10(x,par);
    opt=optimset('Display','off','MaxFunEvals',1e4,'MaxIter',1e4,'TolFun',...
        1e-6,'TolX',1e-6,'Jacobian','on'); %
    [xsol,FVAL,EXITFLAG,OUTPUT,J]=fsolve(funb,y,opt);
%     checkneg2=find(xsol<0);
%     if ~isempty(checkneg2)
%         warning('negative fsolve solution')
%         simdata=feval(@datafile_prv10,tend,y,paramvec);
%     else
%         simdata=feval(@datafile_prv10,tend,xsol,paramvec);
%     end
elseif l==6
    options = odeset('RelTol',1e-3,'AbsTol',atol,'NonNegative', 1:length(ic),...
        'Stats','off','Jacobian',@v10sjac); %,'OutputSel',[21,22],'OutputFcn',@odeplot
    [t,xout] = ode23tb(@prv10,[0 tf],ic,options,par);
    xtemp =xout(end,:);
    %J=0;
    funb=@(x) fsolprv10(x,par);
    opt=optimset('Display','off','MaxFunEvals',1e4,'MaxIter',1e4,'TolFun',...
        1e-6,'TolX',1e-6,'Jacobian','on'); %
    [xsol,FVAL,EXITFLAG,OUTPUT,J]=fsolve(funb,xtemp,opt);
    return
end

function [out flag new_data]=rhs(t,y,data)
par=data.p;
out=prv10(t,y,par);
flag=0;
new_data=[];
return
%Jacobian definition
function [J,flag,new_data]=jv10(t,y,fy,data)
%% State Variables
p=data.p;
J=v10densesjac(t,y,p);
flag=0;
new_data=[];
return
