function [output,output_pre,output_post] = extreme_GHK_elecoPoration_resistance(cell_input)

no_Pairs = size(cell_input,2);

pre_res = [];
post_res = [];
output_pre = [];
output_post = [];
output = [];

for i = 1:no_Pairs
    cap_Pair = cell_input{i};
    [cap1,cap1Num] = LoadExperiments( cap_Pair(1), 1, [0 16] );
    [cap2,cap2Num] = LoadExperiments( cap_Pair(2), 1, [0 16] );
    for j = 1:cap1Num
        pre_res(j,1) = cap1(j).getResistance;
        pre_res(j,2) = cap1(j).getReservoirConc();
    end
    for k = 1:cap2Num
        post_res(k,1) = cap2(k).getResistance;
        post_res(k,2) = cap2(k).getReservoirConc();
    end
    res = post_res./pre_res;
    output_pre = catpad(2,output_pre,pre_res,'padval',0);
    output_post = catpad(2,output_post,post_res,'padval',0);
    output = catpad(2,output,res,'padval',0);
end

end