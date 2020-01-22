kcl_1mM = [1014]; %change descripon
kcl_10mM = [180];
kcl_100mM = [155,156,160:163,165:168,198,206:207,238:239,870];
kcl_1M = [181,873,875,878,883,886,888,890,891,896,900,901,903,908,911,919,923];

[params_1mM,output_1mM] = extreme_GHK_Fit(kcl_1mM,'KCl','cation');
[params_10mM,output_10mM] = extreme_GHK_Fit(kcl_10mM,'KCl','cation');
[params_100mM,output_100mM] = extreme_GHK_Fit(kcl_100mM,'KCl','cation');
[params_1M,output_1M] = extreme_GHK_Fit(kcl_1M,'KCl','cation');