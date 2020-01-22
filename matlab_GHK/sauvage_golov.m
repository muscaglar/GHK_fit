function [s_ghk] = sauvage_golov(conc_1,conc_2,Vnernsts,z_ca,z_an)

F = 96485.33289; % C mol-1
R = 8.3144598; % m2 kg s-2 K-1 mol-1
T = 298; % K
s_ghk = [];
cons = F/(R*T);

for i = 1:length(Vnernsts)
    
    c_high = conc_2(i); %c_res max(conc_2(i),conc_1(i));
    c_low = conc_1(i); %c_cap min(conc_2(i),conc_1(i));
    
    Cin = min(conc_2(i),conc_1(i));%conc_1(i); %c_res
    Cout = max(conc_2(i),conc_1(i));%conc_2(i); %c_cap
    
    Vnernst = Vnernsts(i)*0.001;
    
%    term_1  = (( (z_an(1)^2) * ((c_low * abs(z_an(1)))/abs(z_ca(1))) - ((c_high * abs(z_an(1)))/abs(z_ca(1)))*exp(Vnernst*z_an(1)*cons) )) / (1 - exp(Vnernst*z_an(1)*cons));
%    term_2 = ((z_ca(1)^2) * (1-exp(Vnernst*z_ca(1)*cons))) / (c_low - c_high*exp(Vnernst*z_ca(1)*cons));
%    s_ghk = [s_ghk, -1*term_1 * term_2];

s_ghk_current = (abs(z_an(1)) / z_ca(1)) * (  (1 - exp(-z_ca(1)*Vnernst*cons) ) / (Cin -Cout*exp(-z_ca(1)*Vnernst*cons) ))* ( (Cin -Cout*exp(-z_an(1)*Vnernst*cons)) / ( 1 - exp(-z_an(1)*Vnernst*cons)) );
s_ghk = [s_ghk, s_ghk_current];

end

s_ghk = s_ghk';

end
