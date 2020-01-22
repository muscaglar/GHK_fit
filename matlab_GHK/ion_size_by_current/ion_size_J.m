function [names, results] = ion_size_J()

cecl3_100mM = [];
caps_CeCl3;
hfcl4_100mM = [];
caps_HfCl4;
K3PO4_100mM = [];
caps_K3PO4;
kcl_100mM = [];
caps_KCl;
lacl3_100mM = [];
caps_LaCl3;
licl_100mM = [];
caps_LiCl;
MgCl2_100mM = [];
caps_MgCl2;
zrcl4_100mM = [];
caps_ZrCl4;

catalogue_pH;

z = 1;

results = [];
names = [];

for cap = cecl3_100mM
    ce_cecl3 = 1.196; %ASSUMING, IX and 3 charge.
    cl_cecl3 = 1.81;
    [allResConcs,allCapConcs,total_aVoff,total_aIoff,~, ~,~, ~,~] = ca_selectivity(cap);%plot,custom_min,custom_max);
    concs = size_check(allResConcs,allCapConcs,total_aVoff,total_aIoff);
    conc_index = size(concs,1);
    for j = 1:conc_index
        res_sigma = CeCl3(concs(j,1)*1000);
        res_sigma = res_sigma(2);
        cap_sigma = CeCl3(concs(j,2)*1000);
        cap_sigma = cap_sigma(2);
        results(z,:) = [ce_cecl3,cl_cecl3,res_sigma,cap_sigma,concs(j,:)];
        names{z} = ['cecl3'];
        z = z + 1;
    end
end

for cap = hfcl4_100mM
    hf_hfcl4 = 0.58;%IV
    cl_hfcl4 = 1.81;
    [allResConcs,allCapConcs,total_aVoff,total_aIoff,~, ~,~, ~,~] = ca_selectivity(cap);%plot,custom_min,custom_max);
    concs = size_check(allResConcs,allCapConcs,total_aVoff,total_aIoff);
    conc_index = size(concs,1);
    for j = 1:conc_index
        res_sigma = HfCl4(concs(j,1)*1000);
        res_sigma = res_sigma(2);
        cap_sigma = HfCl4(concs(j,2)*1000);
        cap_sigma = cap_sigma(2);
        results(z,:) = [hf_hfcl4,cl_hfcl4,res_sigma,cap_sigma,concs(j,:)];
        names{z} = ['hfcl4'];
        z = z + 1;
    end
end

for cap = K3PO4_100mM
    k_k3po4 = 1.51; %VIII  ASSUMPTION
    po4_k3po4 = 2.38;
    [allResConcs,allCapConcs,total_aVoff,total_aIoff,~, ~,~, ~,~] = ca_selectivity(cap);%plot,custom_min,custom_max);
    concs = size_check(allResConcs,allCapConcs,total_aVoff,total_aIoff);
    conc_index = size(concs,1);
    for j = 1:conc_index
        res_sigma = K3PO4(concs(j,1)*1000);
        res_sigma = res_sigma(2);
        cap_sigma = K3PO4(concs(j,2)*1000);
        cap_sigma = cap_sigma(2);
        results(z,:) = [k_k3po4,po4_k3po4,res_sigma,cap_sigma,concs(j,:)];
        names{z} = ['K3PO4'];
        z = z + 1;
    end
    
end

for cap = kcl_100mM
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

for cap = lacl3_100mM
    la_lacl3 = 1.1; %VII
    cl_lacl3 = 1.81;
    [allResConcs,allCapConcs,total_aVoff,total_aIoff,~, ~,~, ~,~] = ca_selectivity(cap);%plot,custom_min,custom_max);
    concs = size_check(allResConcs,allCapConcs,total_aVoff,total_aIoff);
    conc_index = size(concs,1);
    for j = 1:conc_index
        res_sigma = LaCl3(concs(j,1)*1000);
        res_sigma = res_sigma(2);
        cap_sigma = LaCl3(concs(j,2)*1000);
        cap_sigma = cap_sigma(2);
        results(z,:) = [la_lacl3,cl_lacl3,res_sigma,cap_sigma,concs(j,:)];
        names{z} = ['LaCl3'];
        z = z + 1;
    end
    
end

for cap = licl_100mM
    li_licl = 0.58;
    cl_licl = 1.81;
    [allResConcs,allCapConcs,total_aVoff,total_aIoff,~, ~,~, ~,~] = ca_selectivity(cap);%plot,custom_min,custom_max);
    concs = size_check(allResConcs,allCapConcs,total_aVoff,total_aIoff);
    conc_index = size(concs,1);
    for j = 1:conc_index
        res_sigma = LiCl(concs(j,1)*1000);
        res_sigma = res_sigma(2);
        cap_sigma = LiCl(concs(j,2)*1000);
        cap_sigma = cap_sigma(2);
        results(z,:) = [li_licl,cl_licl,res_sigma,cap_sigma,concs(j,:)];
        names{z} = ['LiCl'];
        z = z + 1;
    end
    
end
for cap = MgCl2_100mM
    mg_mgcl2 = 0.72; %VI
    cl_mgcl2 = 1.81;
    [allResConcs,allCapConcs,total_aVoff,total_aIoff,~, ~,~, ~,~] = ca_selectivity(cap);%plot,custom_min,custom_max);
    concs = size_check(allResConcs,allCapConcs,total_aVoff,total_aIoff);
    conc_index = size(concs,1);
    for j = 1:conc_index
        res_sigma = MgCl2(concs(j,1)*1000);
        res_sigma = res_sigma(2);
        cap_sigma = MgCl2(concs(j,2)*1000);
        cap_sigma = cap_sigma(2);
        results(z,:) = [mg_mgcl2,cl_mgcl2,res_sigma,cap_sigma,concs(j,:)];
        names{z} = ['MgCl2'];
        z = z + 1;
    end
    
end
for cap = zrcl4_100mM
    zr_zrcl4 = 0.84; %VIII
    cl_zrcl4 = 1.81;
    [allResConcs,allCapConcs,total_aVoff,total_aIoff,~, ~,~, ~,~] = ca_selectivity(cap);%plot,custom_min,custom_max);
    concs = size_check(allResConcs,allCapConcs,total_aVoff,total_aIoff);
    conc_index = size(concs,1);
    for j = 1:conc_index
        res_sigma = ZrCl4(concs(j,1)*1000);
        res_sigma = res_sigma(2);
        cap_sigma = ZrCl4(concs(j,2)*1000);
        cap_sigma = cap_sigma(2);
        results(z,:) = [zr_zrcl4,cl_zrcl4,res_sigma,cap_sigma,concs(j,:)];
        names{z} = ['ZrCl4'];
        z = z + 1;
    end
    
end
names = names';

results(:,9)=results(:,5)./results(:,6);
results(:,10)=results(:,3)./results(:,4);
results(:,11) = results(:,8)./results(:,9);
results(:,12)=results(:,8)./results(:,10);

end