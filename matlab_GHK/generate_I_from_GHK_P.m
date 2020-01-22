%Vn = [-181.0265,-139.94483,-119.943,-99.92863,-79.96458,-60.00256,-39.98341,-19.99254,0.00231,20.00677,39.96925,59.96922,79.97417,99.87571,119.92236,139.88843,179.93154,199.92075];
%Vn = [-119.97717,-80.02826,-39.98455,-0.02987,39.97522,79.94262,119.91446,159.875];
%Vn = [-179.9215,-159.93369,-139.90579,-119.94321,-99.99221,-80.18361,-60.053,-40.25225,-20.09182,-0.05253,19.97999,39.99232,59.98674,79.92632,99.9571,119.94305,139.94018,159.93621,179.90647];
%Vn = [-140.15398,-119.98184,-100.16372,-79.93571,-59.9744,-39.95639,-20.01041,-0.04086,19.98035,39.95611,59.9433,79.94908,99.92858,119.8913,139.93908,159.91577,179.92582];
function [output] = generate_I_from_GHK_P(P,Vnerst,Co,Ci)

Vn = -200:1:200;
Vn = Vn.*-0.001;

max_Vnernst = 59;

%HfCl4

zan = -1;
zcat = 4;

I = [];
I_an = [];
I_cat = [];
I_t = [];

Pcat = 10;
Pan = Pcat*P;

R = 8.3144621;
F = 96485.3329;
T = 300;

Co_an  = Co*(abs(zcat/zan));
Ci_an  = Ci*(abs(zcat/zan));
Co_cat = Co*(abs(zan/zcat));
Ci_cat = Ci*(abs(zan/zcat));

for V = Vn
    an_expTerm = exp(-zan*F*V*(1/R)*(1/T));
    an_I_temp = (Pan*zan*zan*V*F*F)*(Co_an-(Ci_an*an_expTerm))*(1/(R*T*(1-an_expTerm)));
    I_an = [I_an, an_I_temp];
    
    cat_expTerm = exp(-zcat*F*V*(1/R)*(1/T));
    cat_I_temp = (Pcat*zcat*zcat*V*F*F)*(Co_cat-(Ci_cat*cat_expTerm))*(1/(R*T*(1-cat_expTerm)));
    I_cat = [I_cat, cat_I_temp];
end

I_an = I_an'.*10E-9;
I_cat = I_cat'.*10E-9;
Vn = Vn.*(1000);
Vn = Vn.*(Vnerst/max_Vnernst);
output = [ Vn' I_an I_cat];

% I_tot = I_an + I_cat;
% I_stepSize = I_tot(10)-I_tot(11);
% zeroIndex = isalmost(I_tot,0,I_stepSize/4);
% VNernst = mean(Vn(zeroIndex))*1000

I_tot = I_an + I_cat;
zeroIndex = isalmost(Vn,0,0.001);
INernst = nanmean(I_tot(zeroIndex));

end
