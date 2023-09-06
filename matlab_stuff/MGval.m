function MGval = MGval(x,Evib,gamma,data,dP,Kp,q,Vb)
% Mie gruneisen evaluation for K0, Kp and Vo for a given gamma and Evib

T = data(:,2);
P = data(:,1);
V = data(:,3);

if (length(x)==4)
    K0 = x(1);
    Vo = x(4);
    Kp = x(2);
    Kpp = x(3);
    Km = [Kpp, Kp, K0];
elseif (length(x)==3)
    K0 = x(1);
    Vo = x(3);
    Kp = x(2);
    Km = [Kp, K0];
else
    K0 = x(1);
    vo = x(2);
    Km = [Kp, K0];
end

Pthermal=1e-6*gamma*V.^-1.*(V/Vb).^q.*Evib';
Pdc=finite_strain(V.^-1,Vo^-1,Km,[1 1]); % residual data - Cold compression curve
Pd=Pdc+Pthermal; % ideal data
%MGval = sum((P-Pd).^2 ./ (dP .* dP)); % residual data - Cold compression curve
MGval = sum((P-Pd).^2);