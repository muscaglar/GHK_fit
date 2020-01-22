function [log_Conc_Ratios,anion_s_ghk,cation_s_ghk,c_res] = sauvage_ideal_GHK(z_an,z_ca)

%Scenarios:
%Completely anion selective
%Completely cation selective

c_cap = [];
c_res = [];
c = [0.0001, 0.001, 0.01,0.1,1,2,4];
for j = 1:length(c)
    for k = 1:length(c)
        c_cap = [c_cap,c(j)];
    end
    for i = 1:length(c)
        c_res = [c_res,c(i)];
    end
end

c_cap = c_cap';
c_res = c_res';

log_Conc_Ratios = log10(c_cap ./ c_res );

[ anion_Nernst ] = percent_Nernst( log_Conc_Ratios, z_an );
[ cation_Nernst ] = percent_Nernst( log_Conc_Ratios, z_ca );

[anion_s_ghk] = sauvage_golov(c_res,c_cap,anion_Nernst*1000,z_ca,z_an);
[cation_s_ghk] = sauvage_golov(c_res,c_cap,cation_Nernst*1000,z_ca,z_an);

figure;
scatter(log_Conc_Ratios,cation_s_ghk);
hold on;
scatter(log_Conc_Ratios,anion_s_ghk);
set(gca,'xscale','log');
set(gca,'yscale','log');

figure;
scatter(c_res,cation_s_ghk);
hold on;
scatter(c_res,anion_s_ghk);
set(gca,'xscale','log');
set(gca,'yscale','log');

end