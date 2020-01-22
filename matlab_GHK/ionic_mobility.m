function [conc,specific_conductance] = ionic_mobility()

cecl3_100mM = []; caps_CeCl3;
hfcl4_100mM = []; caps_HfCl4;
K3PO4_100mM = []; caps_K3PO4;
kcl_combined = [];caps_KCl;
lacl3_100mM = []; caps_LaCl3;
licl_100mM = [];  caps_LiCl;
MgCl2_100mM = []; caps_MgCl2;
zrcl4_100mM = []; caps_ZrCl4;

catalogue_pH;

reuslts = [];
names = [];
z = 1;

for cap = kcl_combined
    k_kcl = 1.51; %VIII
    cl_kcl = 1.81;
    [allResConcs,allCapConcs,total_aVoff,total_aIoff,~, ~,~, ~,~] = ca_selectivity(cap);%plot,custom_min,custom_max);
    concs = size_check(allResConcs,allCapConcs,total_aVoff,total_aIoff);
    conc_index = size(concs,1);
    for j = 1:conc_index
        res_sigma = KCl(concs(j,1)*1000);
        res_sigma = res_sigma(2);
        cap_sigma = KCl(concs(j,2)*1000);
        cap_sigma = cap_sigma(2);
        results(z,:) = [k_kcl,cl_kcl,res_sigma,cap_sigma,concs(j,:)];
        names{z} = ['KCl'];
        z = z + 1;
    end  
end

index = results(:,5)==results(:,6);


conc = results(index,5);
conc_stoich = conc;
specific_conductance = results(index,3)./100;
specific_conductivity = 1000*specific_conductance ./ conc_stoich;
close all;

scatter(conc_stoich,specific_conductivity);


%CatX AnY
%KCl
%X= 1, Y = 1
 
%specific_conductance = x * lambda_cat*conc_stoich + y*lambda_an*conc_stoich;


%specific_conductivity = x * cat_mobility + y*an_mobility; %solve simultaneously
end