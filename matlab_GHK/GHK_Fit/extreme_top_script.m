caps_KCl_LaCl3;
extreme_GHK_Fit_GO_curve_fit('KCl_LaCl3','Cation',KCl_100mM_LaCl3_10mM,'10mM',KCl_100mM_LaCl3_100mM,'100mM',KCl_100mM_LaCl3_1M,'1M',KCl_100mM_LaCl3_2M,'2M');
clear all;

caps_HfCl4;
extreme_GHK_Fit_GO_curve_fit('HfCl4','Anion',hfcl4_100mM,'100mM',hfcl4_10mM,'10mM',hfcl4_1mM,'1mM',hfcl4_1M,'1M');
clear all;

caps_ZrCl4;
extreme_GHK_Fit_GO_curve_fit('ZrCl4','Anion',zrcl4_100mM,'100mM',zrcl4_10mM,'10mM',zrcl4_1M,'1M');
clear all;

caps_CeCl3;
extreme_GHK_Fit_GO_curve_fit('CeCl3','Anion',cecl3_100mM,'100mM',cecl3_10mM_anion,'10mM_An',cecl3_10mM_both,'10mM',cecl3_1M,'1M',cecl3_2M,'2M',cecl3_3M,'3M');
extreme_GHK_Fit_GO_curve_fit('CeCl3','Cation',cecl3_10mM_cation,'10mM_Ca',cecl3_1mM,'1mM');
clear all;

caps_LaCl3;
extreme_GHK_Fit_GO_curve_fit('LaCl3','Anion',lacl3_100mM,'100mM',lacl3_10mM_anion,'10mM_An',lacl3_1M,'1M',lacl3_2M,'2M',lacl3_3M,'3M',lacl3_4M,'4M');
extreme_GHK_Fit_GO_curve_fit('LaCl3','Cation',lacl3_10mM_cation,'10mM_Ca',lacl3_1mM,'1mM');
extreme_GHK_Fit_GO_curve_fit('ALaCl3','Anion',alacl3_100mM,'100mM',alacl3_1mM,'1mM',alacl3_1M,'1M');
extreme_GHK_Fit_GO_curve_fit('ALaCl3','Cation',alacl3_10mM,'10mM');
clear all;

caps_KCl;
extreme_GHK_Fit_GO_curve_fit('KCl','Cation',kcl_100mM,'100mM',kcl_10mM,'10mM',kcl_1mM,'1mM',kcl_1M,'1M');
clear all;

caps_BKCl
extreme_GHK_Fit_GO_curve_fit('BKCl','Cation',KCl_B_100mM,'100mM',KCl_B_10mM,'10mM',KCl_B_1M,'1M');
clear all;

caps_K3PO4
extreme_GHK_Fit_GO_curve_fit('K3PO4','Cation',K3PO4_100mM,'100mM',K3PO4_10mM,'10mM',K3PO4_1M,'1M');
clear all;

caps_KCl_HfCl4
extreme_GHK_Fit_GO_curve_fit('KCl_HfCl4','Cation',KCl_100mM_HfCl4_1mM,'100mM',KCl_100mM_HfCl4_0_1mM,'100mM',KCl_100mM_HfCl4_0_01mM,'100mM');
clear all;

% caps_KCl_LaCl3;
% extreme_GHK_Fit_GO_test('KCl_LaCl3','Cation',KCl_100mM_LaCl3_10mM,'10mM',KCl_100mM_LaCl3_100mM,'100mM',KCl_100mM_LaCl3_1M,'1M',KCl_100mM_LaCl3_2M,'2M');
% clear all;
% 
% caps_HfCl4;
% extreme_GHK_Fit_GO_test('HfCl4','Anion',hfcl4_100mM,'100mM',hfcl4_10mM,'10mM',hfcl4_1mM,'1mM',hfcl4_1M,'1M');
% clear all;
% 
% caps_ZrCl4;
% extreme_GHK_Fit_GO_test('ZrCl4','Anion',zrcl4_100mM,'100mM',zrcl4_10mM,'10mM',zrcl4_1M,'1M');
% clear all;
% 
% caps_CeCl3;
% extreme_GHK_Fit_GO_test('CeCl3','Anion',cecl3_100mM,'100mM',cecl3_10mM_anion,'10mM_An',cecl3_10mM_both,'10mM',cecl3_1M,'1M',cecl3_2M,'2M',cecl3_3M,'3M');
% extreme_GHK_Fit_GO_test('CeCl3','Cation',cecl3_10mM_cation,'10mM_Ca',cecl3_1mM,'1mM');
% clear all;
% 
% caps_LaCl3;
% extreme_GHK_Fit_GO_test('LaCl3','Anion',lacl3_100mM,'100mM',lacl3_10mM_anion,'10mM_An',lacl3_1M,'1M',lacl3_2M,'2M',lacl3_3M,'3M',lacl3_4M,'4M');
% extreme_GHK_Fit_GO_test('LaCl3','Cation',lacl3_10mM_cation,'10mM_Ca',lacl3_1mM,'1mM');
% extreme_GHK_Fit_GO_test('ALaCl3','Anion',alacl3_100mM,'100mM',alacl3_1mM,'1mM',alacl3_1M,'1M');
% extreme_GHK_Fit_GO_test('ALaCl3','Cation',alacl3_10mM,'10mM');
% clear all;
% 
% caps_KCl;
% extreme_GHK_Fit_GO_test('KCl','Cation',kcl_100mM,'100mM',kcl_10mM,'10mM',kcl_1mM,'1mM',kcl_1M,'1M');
% clear all;
% 
% caps_BKCl
% extreme_GHK_Fit_GO_test('BKCl','Cation',KCl_B_100mM,'100mM',KCl_B_10mM,'10mM',KCl_B_1M,'1M');
% clear all;
% 
% caps_K3PO4
% extreme_GHK_Fit_GO_test('K3PO4','Cation',K3PO4_100mM,'100mM',K3PO4_10mM,'10mM',K3PO4_1M,'1M');
% clear all;
% 
% caps_KCl_HfCl4
% extreme_GHK_Fit_GO_test('KCl_HfCl4','Cation',KCl_100mM_HfCl4_1mM,'100mM',KCl_100mM_HfCl4_0_1mM,'100mM',KCl_100mM_HfCl4_0_01mM,'100mM');
% clear all;


% caps_HfCl4;
% extreme_GHK_Fit_GO('HfCl4',hfcl4_100mM,'100mM',hfcl4_10mM,'10mM',hfcl4_1mM,'1mM',hfcl4_1M,'1M');
% clear all;
% 
% caps_ZrCl4;
% extreme_GHK_Fit_GO('ZrCl4',zrcl4_100mM,'100mM',zrcl4_10mM,'10mM',zrcl4_1M,'1M');
% clear all;

% caps_CeCl3;
% extreme_GHK_Fit_GO('CeCl3',cecl3_100mM,'100mM',cecl3_10mM_anion,'10mM_An',cecl3_10mM_cation,'10mM_Ca',cecl3_10mM_both,'10mM',cecl3_1mM,'1mM',cecl3_1M,'1M',cecl3_2M,'2M',cecl3_3M,'3M');
% clear all;
% 
% caps_LaCl3;
% extreme_GHK_Fit_GO('LaCl3',lacl3_100mM,'100mM',lacl3_10mM_anion,'10mM_An',lacl3_10mM_cation,'10mM_Ca',lacl3_1mM,'1mM',lacl3_1M,'1M',lacl3_2M,'2M',lacl3_3M,'3M',lacl3_4M,'4M');
% extreme_GHK_Fit_GO('ALaCl3',alacl3_100mM,'100mM',alacl3_10mM,'10mM',alacl3_1mM,'1mM',alacl3_1M,'1M');
% clear all;
% 
% caps_KCl;
% extreme_GHK_Fit_GO('KCl',kcl_100mM,'100mM',kcl_10mM,'10mM',kcl_1mM,'1mM',kcl_1M,'1M');
% clear all;
% 
% caps_BKCl
% extreme_GHK_Fit_GO('BKCl',KCl_B_100mM,'100mM',KCl_B_10mM,'10mM',KCl_B_1M,'1M');
% clear all;
% 
% caps_K3PO4
% extreme_GHK_Fit_GO('K3PO4',K3PO4_100mM,'100mM',K3PO4_10mM,'10mM',K3PO4_1M,'1M');
% clear all;