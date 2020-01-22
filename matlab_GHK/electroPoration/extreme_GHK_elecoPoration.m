function [output] = extreme_GHK_elecoPoration(cell_input,saltName,polarity,Origin)

no_Pairs = size(cell_input,2);
output = zeros(no_Pairs,12);

output1 = zeros(no_Pairs,7);
output2 = zeros(no_Pairs,7);
output3 = zeros(no_Pairs,7);
output4 = zeros(no_Pairs,7);

resistances_t = [];

for i = 1:no_Pairs
    cap_Pair = cell_input{i};
    [resistances,pore_sizes,~] = amoureux_selectivityPoreAnalysis(cap_Pair);
    resistances_t =    [resistances_t; resistances];
    [pre_params,~,pre_params_Alt] = extreme_GHK_Fit(cap_Pair(1),saltName,polarity,[],[],[],0);%dirName,plotName,sheetName,Origin)
    [post_params,~,post_params_Alt] = extreme_GHK_Fit(cap_Pair(2),saltName,polarity,[],[],[],0);%dirName,plotName,sheetName,Origin)
    
    close all;
    output1(i,:) = pre_params;
    output2(i,:) = pre_params_Alt;
    
    output3(i,:) = post_params;
    output4(i,:) = post_params_Alt;
    
    output(i,:) = [pore_sizes(3),resistances(3),resistances(6),resistances(9),pre_params(2),pre_params(4),post_params(2),post_params(4),pre_params_Alt(2),pre_params_Alt(4),post_params_Alt(2),post_params_Alt(4)];
end
if(Origin)
    extreme_GHK_Params_Origin(output1,['pre' saltName 'P'],[]);
    extreme_GHK_Params_Origin(output2,['pre' saltName 'PA'],[]);
    extreme_GHK_Params_Origin(output3,['post' saltName 'P'],[]);
    extreme_GHK_Params_Origin(output4,['post' saltName 'PA'],[]);
    extreme_GHK_Params_Origin(output(:,[2:4]),['resistances_' saltName],[]);
end
end


%Need to figure out where these resistances are being calculated from and
%if we can do better by checking Neroli.