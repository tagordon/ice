function sp=sp_F_fit(data,options)
%function to generate a LBF of Helmholz energy using an arbitrary strain metric
% Usage:  sp=sp_F_fit(data,options)
%   where the structure "data" contains:
%                         PV an array of pressure - volume points
%                         VKP an optional array of bulk moduli data
%                         Vo - the volume associated with zero strain
%
%        the structure "options" contains:
%               k defines order (if not supplied, default to 6)
%               Vflg flg ==1 fo input in volumes ==0 for input in inverse V
%               strainflg character string for various strain
%                                  metrics to use.  Must be set, seee
%                                  function getStrain for choices
%               K_weight optional weighting of bulk moduli data (default is no K data)
%               knt:  contains the knt sequence to use in volume or
%                             density units
%               kntflg if set to 1, directly use the knots in options.knt
%               Reg if missing, no regularization, scalar number defines a uniform
%                            number of regularization points between first
%                            and last knot, or a vector defines a set of point location 
%               drv  derivative of energy to minimize in regulation
%               lam  weight of regularization
%    units: currently using GPa and gm/cc - this should be reviewed
%   JMB 2019

% it does not matter whether the independent parameter is density or volume
%  but we need to be consistent. internally use V 

% Preliminaries
% default spline order is 6
if(isfield(options,'k'))
    k=options.k;
else
    k=6;
end

% whether density or volume data are input has to be specified
if(~isfield(options,'Vflg'))
     error('Please set options.Vflg (1 for volume input, 0 for density input')
else  % input specified - act on it
    if options.Vflg==1  % volume data entered
      Vo=data.Vo;
      V=data.PV(:,2);
      knt=options.knt;
    elseif options.Vflg==0  % density data entered
      Vo=data.Vo.^-1;
      V=data.PV(:,2).^-1;
      knt=options.knt.^-1;
    else
        error('Unknown flag for volume input')
    end   
end

% it is annoying but the function collocate requires strains to increase
% thus data and knots (here in volumes) are sorted into reverse order. 
[V,id]=sort(V,'descend'); %order V so that strain increases
P=data.PV(id,1);
knt=sort(knt,'descend'); %order V so that strain increases

% determine if bulk modulus data are to be fit
Kflg=0;
if(isfield(options,'K_weight'))
    Kflg=1;
    if(~isfield(data,'K')), error('Fitting K requires data.K'),end   
    lamK=options.K_weight;
end

% everything below should be consistently in volume units
if Kflg  % only need to sort K if K data are input
    VK=data.K(:,1);
    if options.Vflg==0 
        VK=VK.^-1;
    end
    [VK,id]=sort(VK,'descend'); %order V so that strain increases
    K=data.K(id,2);
end

if (isfield(options,'strainflg'))
    out=getStrains(knt,Vo,options.strainflg); % turn knots in V to knots in strain
    knt=out.f;
else
    error('strain metric needs to be specified in options.strainflg')
end

% an option to hand specify knots in appropriate strain metrics
if(isfield(options,'kntflg'))
    if options.kntflg==1
        knt=options.knt;
    end
end

% determine the strains for P-V and K-V  input
out=getStrains(V,Vo,options.strainflg);
fd=out.f;
fvd=out.fv;

if Kflg
    out=getStrains(VK,Vo,options.strainflg);
    fK=out.f;
    fvK=out.fv;
    f2vK=out.f2v;
end

% add k-fold multiplicity to the first and last knot
knt=[repmat(knt(1),1,k) knt(2:end-1) repmat(knt(end),1,k)];

%*********************************************************************************
% next section deals with regularization (which is not always required).
% Points to regularize are input in options.Reg in either volume or
% density units (as specified by the Vflg option) or as a scalar number 
% of points between the first and last knots 
Regflg=0;
if(isfield(options,'Reg'))
    if length(options.Reg)==1
       nReg=options.Reg;
       fr=linspace(knt(1),knt(end),nReg);
    else
        nReg=length(options.Reg);
        if (isfield(options,'strainflg'))
            if options.Vflg
               out=getStrains(options.Reg,Vo,options.strainflg);
               fr=out.f;
            else
               out=getStrains(options.Reg.^-1,Vo,options.strainflg);
               fr=out.f;
            end
        else
            if options.Vflg
                fr=.5*((options.Reg/Vo).^(-2/3)-1);  
            else
                fr=.5*((options.Reg.^-1/Vo).^(-2/3)-1);   
            end
            fr=sort(fr);
        end
    Rcol=collocate(knt,k,fr,k-1);
    Regflg=1;
    
    lam=options.lam;  % the weight for regularization is set
    drv=1+options.drv; % the derivative to minimize is set. (
    wtReg=norm(Rcol(:,:,drv),1)^-1; % the "1 norm" of the basis functions provides a plausible weight
    end
end

%***********************************************************************************
% get basis functions for data (pressure and bulk moduli:

Pcol=collocate(knt,k,fd,2);

if Kflg
    Kcol=collocate(knt,k,fK,2);
    [~,n,~]=size(Kcol);
    %the linear combination of derivatives that gives K/V:
    C=-Kcol(:,:,3).*(fvK(:).^2*ones(1,n))-Kcol(:,:,2).*(f2vK(:)*ones(1,n)); 
end

% get the basis function for Helmholz energy at Vo or zero strain (need
% energy data since we are "collocating" from hiher order derivatives of
% energy
if (isfield(options,'strainflg'))
    
    if(strcmp(options.strainflg(1:3),'vol'))        
        Fcol=collocate(knt,k,Vo,0);
        % annoying - we previously ordered data for strain metrics that 
        %  increase with pressure while volume decreases and thus the 
        %    pressures need to be flipped (V is flipped by the getstrains
        %    function
         P=P(end:-1:1); 
         V=V(end:-1:1); 
         K=K(end:-1:1); 
    else
         Fcol=collocate(knt,k,0,0);
    end
else       
     Fcol=collocate(knt,k,0,0);   % define zero energy at zero strain - can add actual energy back later
end

%****************************************************************************************
% assemble the basis function and data arrays for the least square fit
% Basis functions in array , data in vector  B

datP=-P./fvd;  % to have basis functions for derivatives of Helmholz energy on the left side
%                  and data on the right P data are "normalized" by the volume derivative of the 
%                  strain metric

% the easy way to normalize the combination of different data that span different orders 
%of magnitude is to divide data by the standard deviation of the data
wtP=1/std(datP); % the weight for pressure data 
if Kflg % data fit of P and K
    if(strcmp(options.strainflg(1:3),'vol'))
        VK=VK(end:-1:1);
    end
    datK=-K./VK; %  The separation of basis functions from "data"
    wtK=1/std(datK); % the normalization weight for K
    %  normalized and weighted basis functions
    A=[Fcol;wtP*Pcol(:,:,2);wtK*lamK*C];
    %  normalized and weighted data
    B=[0;wtP*datP;wtK*lamK*datK];
else % data fit of P only
    A=[Fcol;wtP*Pcol(:,:,2)];
    B=[0;-wtP*P./fvd];
end

if Regflg
    A=[A;lam*wtReg*Rcol(:,:,drv)];
    B=[B;zeros(nReg,1)];
end

%***********************************************************************************

% SOLUTION TIME!
%here it is, the least square solution:
m=(A'*A)\(A'*B);


% *************************************************************************************
% all that is left is assembling results in an output structure
sp.form='B-';
sp.knots=knt;
sp.coefs=m(:)';
sp.number=length(m);
sp.order=k;
sp.dim=1;
sp.label='Helmholz representation in a chosen strain metric';
sp.strainflg=options.strainflg;
sp.coefs_units='Energy unit depends on input: GPa/(g/cm^3) = MJ/kg units';
sp.timestamp=datetime;
sp.rho0=1/Vo;

pc=  -(Pcol(:,:,2)*m(:)).*fvd(:);
if (isfield(options,'strainflg'))
end
sp.Data.rms=sqrt(sum(((pc(:)-P(:))).^2)/length(fd));
sp.Data.rms_n=sqrt(sum(((pc(:)-P(:))./P(:)).^2)/length(fd));
sp.Data.PV=[P(:) V(:) pc(:) ];
if Regflg
damp=Rcol(:,:,drv)*m(:);
sp.Data.damp=sqrt(sum((damp).^2)/length(damp));
end

if Kflg 
    kc=-VK.*(C*m(:));
    nc=length(kc);
    sp.Data.rms=[sp.Data.rms sqrt(sum((kc(:)-K(:)).^2)/nc)];
    sp.Data.K=[VK(:) K(:) kc(:)];
end



    



