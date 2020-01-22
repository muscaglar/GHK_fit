function [] = extreme_GHK_Params_Origin_test(input, sheetName,folderDirectory)

ORG = Matlab2OriginPlot();

ORG.cd_TopLevel();
ORG.cd(folderDirectory);

sheetName = ORG.MatrixToOrigin(input, sheetName);

%1: Y params = [aCapConcs(1),fitresult.C,fitresult.A,gof.sse,gof.adjrsquare,gof.rmse];
%2: D/C
%3: Y Error
%4: X
ORG.ExecuteLabTalk('wks.col1.lname$ = Capillary Concentration;');
ORG.ExecuteLabTalk('wks.col1.type = 4;'); 

ORG.ExecuteLabTalk('wks.col2.lname$ = P Cat;');
ORG.ExecuteLabTalk('wks.col2.type = 2;');

ORG.ExecuteLabTalk('wks.col3.lname$ = P An;');
ORG.ExecuteLabTalk('wks.col3.type = 2;');

ORG.ExecuteLabTalk('wks.col4.lname$ = SSE;');
ORG.ExecuteLabTalk('wks.col4.type = 2;');

ORG.ExecuteLabTalk('wks.col5.lname$ = AdjRSquare;');
ORG.ExecuteLabTalk('wks.col5.type = 2;');

ORG.ExecuteLabTalk('wks.col6.lname$ = RMSE;');
ORG.ExecuteLabTalk('wks.col6.type = 2;');


ORG.Disconnect;

end
% [Cap]
% [Res]
% [Cap]/[Res]
% GHK