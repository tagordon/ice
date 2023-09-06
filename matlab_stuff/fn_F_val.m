function out=fn_F_val(sp,rho)
%function to provide EOS quantities from a spline representation of
%Helmholz energy as a function of a chosen strain metric.  In addition to the
%standard spline quantities, sp needs (1) rho0 and (2) a strain metric added to the structure
% Usage: out=fn_F_val(sp,rho)
% where     rho is density or inverse of whatever volume unit was used in
%               construction of the representation
%             out contains P, K and Kp=dKdP in the pressure units of the
%             construction
% JMB 2019
m=sp.coefs(:);
k=sp.order;
V=rho(:).^-1;
Vo=sp.rho0.^-1;

% determine strain and its derivatives wrt V
if (isfield(sp,'strainflg'))
    out=getStrains(V,Vo,sp.strainflg);
    f=out.f;
    fv=out.fv;
    f2v=out.f2v;
    f3v=out.f3v;
else
   error('The strain metric "sp.strainflg" must be specified')
end

% get Helmholz energy basis functions (that are functions of strain)
% The following gets the function and derivatives of the function spline
% basis functions at data site locations

Fcol=collocate(sp.knots,k,f,k-1);

% Determine energy and  first three derivatives of energy wrt strain
F=Fcol(:,:,1)*m;
F1=Fcol(:,:,2)*m;
F2=Fcol(:,:,3)*m;
if k>3
 F3=Fcol(:,:,4)*m;
else
    F3=0;
end

% in the case of using volumes directly the order of the F's needs to be
% reversed
if (isfield(sp,'strainflg'))
 if (strcmp(sp.strainflg(1:3),'vol'))
      F1=F1(end:-1:1);
      F2=F2(end:-1:1);
      F3=F3(end:-1:1); 
 end
end

% determine pressure and first two derivatives of P wrt volume dP/dV and d2P/dV2 
P=-fv.*F1(:);
P1=-F2(:).*fv.^2-F1(:).*f2v;  %derivative of previous line
P2=-F3(:).*fv.^3 -2*F2.*fv.*f2v -F2.*fv.*f2v - F1.*f3v; %derivative of previous line

% calculate energy, pressure, bulk modulus and its first derivative
out.P=P;
out.F=F;
out.K=-V.*P1;
out.Kp=-(V.*P2./P1+1);


