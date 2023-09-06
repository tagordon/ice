function sp_out=NL_F_Fitter(sp,Data,options)
% Function for non-linear optimization of dK/dP using P-V and K-V data as
% well
% Usage:   sp_out=NL_F_Fitter(sp,Data,options)
%    where a starting model is "sp" and data are in structure "Data" and
%    fitting options are in "options"
% see sp_F_fitter for examples of data and options.
% JMB 2019

% put important basis functions and data in structure "Cntrl"
Cntrl=mk_EOSCntrl(sp,Data,options);

sp_out=sp;

mo=sp.coefs;
i=1;
devmax=1;
% from experience 5 or so steps seems to converge. 
while (i<10 & devmax>1e-4)
    i=i+1;
% Here is the one-stop, single step non-linear optimizer.  
    sp_out=NL_step(sp_out,Cntrl);
    mn=sp_out.coefs;
    devmax=max(abs(mn-mo));
    mo=mn;
end



function sp_out=NL_step(sp,Cntrl)

% for nonlinear fitting we need to calculate derivatives wrt to changes in
% the model parameters.  Calculate base model and base Regularization:
mdl0=mdl_EOS(sp,Cntrl);
Reg1=Cntrl.BF.Reg*sp.coefs(:);

dF=1e-10; % set increment of  energy for derivative calculations.

% set up matrixes for the derivatives of the data and regularization
% wrt model parameters
%[np,~]=size(Cntrl.BF.P1);
Adat=zeros(length(mdl0.devs),sp.number);
Areg=zeros(length(Cntrl.BF.Reg(:,1)),sp.number);


% OK - step through and perturb every model parameter 
% and form finite difference derivatives 
for i=1:sp.number
      sptmp=sp;
      sptmp.coefs(i)=sptmp.coefs(i)+dF; % perturb a model parameter
      mdl1=mdl_EOS(sptmp,Cntrl);
% derivatives of all data wrt to each model parameter
      Adat(:,i)=Cntrl.data.weight.*(mdl1.y-mdl0.y)/dF; 
%Derivative of regularization wrt each model parameter:
       Reg2=Cntrl.BF.Reg*sptmp.coefs(:);
       Areg(:,i)   = (Reg2(:)  - Reg1(:))/dF;
end

% assemble the normal equations
A=[    Adat
        Cntrl.BF.Regwt*Cntrl.BF.Reglam*Areg
       ]; 

 B = [mdl0.devs(:).*Cntrl.data.weight
       -Cntrl.BF.Regwt*Cntrl.BF.Reglam*Reg1(:)
       ];
   
%reality check
if not(isempty(find(isnan(B), 1)))
    error('You have NaNs in the B matrix')
end
if not(isempty(find(isnan(A), 1)))
    error('You have NaNs in the A matrix')
end

% and invert for the changes in the model parameters
 dmod = (A'*A)\(A'*B); 
     sp_out=sp;
     sp_out.coefs=sp.coefs+dmod(:)';

mdl_out=mdl_EOS(sp_out,Cntrl);
mdl_out.rms

function mdl=mdl_EOS(sp,Cntrl)
% Function to return model predictions for a b spline "sp" with data sites
% as described in "Cntrl".  Usage:
%  yc=EOS(sp,DataStrc)
% where yc is a vector of predictions in the order
%   [  F(Vo)
%      -pressures/f
%      K/V (if provided
%      Kp at P essentially at infinity
%      ];
% units notes: sp has knots in strain 


% the structure CntrlStrct is created by
% Cntrl.BF.Fcol
% sp.coefs(:)
F=Cntrl.BF.Fcol*sp.coefs(:);
Povf=Cntrl.BF.Pcol*sp.coefs(:);
KovV=Cntrl.BF.Kcol*sp.coefs(:);


P1 = Cntrl.BF.P1*sp.coefs(:);
P2 =Cntrl.BF.P2 *sp.coefs(:);

% calculate energy, pressure, bulk modulus and its first derivative

Kp=-(Cntrl.data.Vkp.*P2./P1+1);

mdl.F=F;
mdl.Povf=Povf;
mdl.KovV=KovV;
mdl.Kp=Kp;

mdl.y=[F;Cntrl.data.wtP*Povf;Cntrl.data.wtK*KovV;Cntrl.data.wtkp*Kp];
mdl.devs=-[F;Cntrl.data.wtP*(Povf-Cntrl.data.datP);Cntrl.data.wtK*(KovV-Cntrl.data.datK);Cntrl.data.wtkp*(Kp-Cntrl.data.datKp)];
mdl.rms=sqrt(sum(mdl.devs.^2/length(mdl.devs)));




function Cntrl=mk_EOSCntrl(sp,Data,options)

%needed: prepare: Fdat, datP wtP  datK  wtK and datKp wtKp Kpcol Rcol wtR
% and pieces to calculate model

knt=sp.knots;
k=sp.order;

Pd=Data.PV(end:-1:1,1);
Vd=Data.PV(end:-1:1,2);
Vo=Data.Vo;

Vkp=Data.Kp(:,1);
Kp=Data.Kp(:,2);

Kflg=0;
if(isfield(options,'K_weight'))
    Kflg=1;
   if(~isfield(Data,'K')), error('Fitting K requires data.K'),end   
   lamK=options.K_weight;
   Vk=Data.K(end:-1:1,4);
   K=Data.K(end:-1:1,2);
end

if (isfield(sp,'strainflg'))
    outd=getStrains(Vd,Vo,sp.strainflg); % turn knots in V to knots in strain
    outo=getStrains(Vo,Vo,sp.strainflg);
    fd=outd.f;
    Pcol=collocate(knt,k,fd,1);
    Pcol=Pcol(:,:,2);
    Fcol=collocate(knt,k,outo.f,0);
    wtP=norm(Pcol,1)^-1;
    if Kflg
        out=getStrains(Vk,Vo,sp.strainflg);
        fk=out.f;
        fvk=out.fv;
        f2vk=out.f2v;
        f3vk=out.f3v;
        Kcol=collocate(knt,k,fk,3);
        [~,n,~]=size(Kcol);
    %the linear combination of derivatives that gives K/V:
         C=-Kcol(:,:,3).*(fvk(:).^2*ones(1,n))-Kcol(:,:,2).*(f2vk(:)*ones(1,n)); 
    end    
else
    error('strain metric needs to be specified in sp.strainflg')
end


regflg=0;
if(isfield(options,'Reg'))
    regflg=1;
     out=getStrains(options.Reg.^-1,Vo,sp.strainflg);
     fr=out.f;
     lam=options.lam;  
     Rcol=collocate(knt,k,fr,k-1);    
     drv=1+options.drv;
     Rcol=Rcol(:,:,drv);
end

datP=-Pd./outd.fv;
wtP=1/std(datP); % the weight for pressure data 

if Kflg % data fit of P and K
    if(strcmp(sp.strainflg(1:3),'vol'))
        Vk=Vk(end:-1:1);
        K=K(end:-1:1);
    end
    datK=-K./Vk; %  The separation of basis functions from "data"
    wtK=1/std(datK);
end

out=getStrains(Vkp,Vo,options.strainflg);

f=out.f;
fv=out.fv(:);
f2v=out.f2v(:);
f3v=out.f3v(:);
Kpcol=collocate(knt,k,out.f,3);
[~,n,~]=size(Kpcol);
P1=-Kpcol(:,:,3).*(fv.^2*ones(1,n))-Kpcol(:,:,2).*(f2v*ones(1,n)); 
P2=-Kpcol(:,:,4).*(fv.^3*ones(1,n)) -2*Kpcol(:,:,3).*((fv.*f2v)*ones(1,n)) -Kpcol(:,:,3).*((fv.*f2v)*ones(1,n)) - Kpcol(:,:,2).*(f3v*ones(1,n));


Cntrl.data.datP=datP;
Cntrl.data.Vd=Vd;
Cntrl.data.wtP=wtP;
nP=length(datP);

Cntrl.data.datK=datK;
Cntrl.data.Vk=Vk;
Cntrl.data.wtK=wtK;
nK=length(K);

Cntrl.data.datKp=Kp;
Cntrl.data.Vkp=Vkp;
Cntrl.data.wtkp=1/mean(Kp);
nKp=length(Kp);

if Kflg
Cntrl.data.weight=[1;ones(nP,1);options.K_weight*ones(nK,1);options.Kp_weight*ones(nKp,1)];
else
  Cntlr.data.weight=[1;ones(nP,1);options.Kp_weight*ones(nKp,1)];  
end


Cntrl.BF.Fcol=Fcol;
Cntrl.BF.Pcol=Pcol;
Pwt=1/norm(Pcol,1);
Cntrl.BF.Kcol=C;
Kwt=1/norm(C,1);
Cntrl.BF.P1=P1;
Cntrl.BF.P1wt=1/norm(P1,1);
Cntrl.BF.P2=P2;
Cntrl.BF.P2wt=1/norm(P2,1);

Cntrl.BF.Reg=Rcol;
Cntrl.BF.Regwt=1/norm(Rcol,1);
Cntrl.BF.Reglam=options.lam;

