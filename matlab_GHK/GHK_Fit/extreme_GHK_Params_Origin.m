function [] = extreme_GHK_Params_Origin(input, sheetName,folderDirectory)

ORG = Matlab2OriginPlot();

ORG.cd_TopLevel();
ORG.cd(folderDirectory);

sheetName = ORG.MatrixToOrigin(input, sheetName);

%1: Y
%2: D/C
%3: Y Error
%4: X


ORG.ExecuteLabTalk('wks.col1.lname$ = Capillary Concetration (M);');
ORG.ExecuteLabTalk('wks.col1.type = 4;'); 

ORG.ExecuteLabTalk('wks.col2.lname$ = P Cat;');
ORG.ExecuteLabTalk('wks.col2.type = 1;');

ORG.ExecuteLabTalk('wks.col3.lname$ = P Cat Err;');
ORG.ExecuteLabTalk('wks.col3.type = 3;');

ORG.ExecuteLabTalk('wks.col4.lname$ = P An;');
ORG.ExecuteLabTalk('wks.col4.type = 1;');

ORG.ExecuteLabTalk('wks.col5.lname$ = P An Err;');
ORG.ExecuteLabTalk('wks.col5.type = 3;');

ORG.ExecuteLabTalk('wks.col6.lname$ = P cat / P an;');
ORG.ExecuteLabTalk('wks.col6.type = 1;');

ORG.ExecuteLabTalk('wks.col7.lname$ = P an / P cat;');
ORG.ExecuteLabTalk('wks.col7.type = 1;');

ORG.Disconnect;

end
% [Cap]
% [Res]
% [Cap]/[Res]
% GHK