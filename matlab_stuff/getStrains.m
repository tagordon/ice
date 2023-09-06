function out=getStrains(V,Vo,strainflg)
%

% note that strain measures increase with pressure while V decreases.
% Thus, here the volume vector is reversed for the case "volume" - this
% requires attention when used.
VovVo=V/Vo;
Vo_ovV=Vo*V.^-1;

switch strainflg
    case 'Eulerian'
            f=1/2*(VovVo).^(-2/3) - 1/2;
            fv=-1/(3*Vo)*(VovVo).^(-5/3);
            f2v=5/(9*Vo^2)*(VovVo).^(-8/3);
            f3v=-40/(27*Vo^3)*(VovVo).^(-11/3);
    case 'log'
            f=-1/3*log(VovVo);
            fv=-1/3*V.^-1;
            f2v=1/3*V.^-2;
            f3v=-2/3*V.^-3;
     case 'volume'
         f=V(end:-1:1);
         fv=1;
         f2v=0;
         f3v=0;
    case 'comp'
        f=-(V - Vo)./V;
        fv=(V - Vo)./V.^2 - V.^-1;
        f2v=2*V.^-2 - (2*(V - Vo)).*V.^-3;
        f3v=(6*(V - Vo)).*V.^-4 - 6*V.^-3;
    case 'Vinet' % cube root of Vo/V
        f=(Vo_ovV).^(1/3)-1;
        fv=-1/(3*Vo)*(Vo_ovV).^(4/3);
        f2v=4/(9*Vo^2)*(Vo_ovV).^(7/3);
        f3v=-28/(27*Vo^3)*(Vo_ovV).^(10/3);
    otherwise
       error('strain uknknown')
      
        
end
out.f=f;
out.fv=fv;
out.f2v=f2v;
out.f3v=f3v;

            
