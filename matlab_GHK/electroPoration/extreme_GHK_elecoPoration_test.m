function [output] = extreme_GHK_elecoPoration_test(cell_input,saltName,polarity,Origin)

no_Pairs = size(cell_input,2);
output = zeros(no_Pairs,30);


resistances_t = [];

for i = 1:no_Pairs
    cap_Pair = cell_input{i};
    [resistances,pore_sizes,~] = amoureux_selectivityPoreAnalysis(cap_Pair);
    resistances_t =    [resistances_t; resistances];
    [pre_params,~] = extreme_GHK_Fit_test(cap_Pair(1),saltName,polarity,num2str(i),num2str(i),num2str(i),0);%dirName,plotName,sheetName,Origin)
    [post_params,~] = extreme_GHK_Fit_test(cap_Pair(2),saltName,polarity,num2str(i),num2str(i),num2str(i),0);%dirName,plotName,sheetName,Origin)
    
    close all;
    
    output(i,:) = [pore_sizes(3),resistances(3),resistances(6),resistances(9),pre_params,post_params];
end

if(Origin)
    extreme_GHK_Params_Origin_test(output,['pre' saltName 'P'],[]);
end

end


%Need to figure out where these resistances are being calculated from and
%if we can do better by checking Neroli.