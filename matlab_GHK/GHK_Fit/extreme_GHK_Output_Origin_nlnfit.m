function [] = extreme_GHK_Output_Origin_nlnfit(input,sheetName,folderDirectory)

ORG = Matlab2OriginPlot();

ORG.cd_TopLevel();
ORG.cd(folderDirectory);


sheetName = ORG.MatrixToOrigin(input, sheetName);

%1: Y
%2: D/C
%3: Y Error
%4: X

ORG.ExecuteLabTalk('wks.col1.lname$ = Capillary Concentration;');
ORG.ExecuteLabTalk('wks.col1.type = 4;'); 

ORG.ExecuteLabTalk('wks.col2.lname$ = P Cat;');
ORG.ExecuteLabTalk('wks.col2.type = 1;'); 

ORG.ExecuteLabTalk('wks.col3.lname$ = P Cat Err;');
ORG.ExecuteLabTalk('wks.col3.type = 3;');

ORG.ExecuteLabTalk('wks.col4.lname$ = P Cat pValue;');
ORG.ExecuteLabTalk('wks.col4.type = 2;');

ORG.ExecuteLabTalk('wks.col5.lname$ = P Cat tStat;');
ORG.ExecuteLabTalk('wks.col5.type = 2;');

ORG.ExecuteLabTalk('wks.col6.lname$ = P An;');
ORG.ExecuteLabTalk('wks.col6.type = 1;');

ORG.ExecuteLabTalk('wks.col7.lname$ = P An Err;');
ORG.ExecuteLabTalk('wks.col7.type = 3;');

ORG.ExecuteLabTalk('wks.col8.lname$ = P An pValue;');
ORG.ExecuteLabTalk('wks.col8.type = 2;');

ORG.ExecuteLabTalk('wks.col9.lname$ = P An tStat;');
ORG.ExecuteLabTalk('wks.col9.type = 2;');

ORG.ExecuteLabTalk('wks.col10.lname$ = Pcat/Pan;');
ORG.ExecuteLabTalk('wks.col10.type = 2;');

ORG.ExecuteLabTalk('wks.col11.lname$ = Pcat/Pan Error;');
ORG.ExecuteLabTalk('wks.col11.type = 3;');

ORG.ExecuteLabTalk('wks.col12.lname$ = Pan/Pcat;');
ORG.ExecuteLabTalk('wks.col12.type = 1;');

ORG.ExecuteLabTalk('wks.col13.lname$ = Pan/Pcat Error;');
ORG.ExecuteLabTalk('wks.col13.type = 3;');

ORG.CreateGraphPage('PcatPan');

ORG.ExecuteLabTalk(['plotxy iy:=[' sheetName ']Sheet1!(1,10,11) plot:=202 ogl:=[' ORG.ActiveGraphName ']' num2str(ORG.CurrentLayer) '!;'] );
colour = 'green';
ORG.CurrentPlotNo = ORG.NoPlots + 1;
ORG.NoPlots = ORG.NoPlots + 1;
ORG.ActivatePage(ORG.ActiveGraphName);
ORG.ExecuteLabTalk(['layer.plot = ' num2str(ORG.CurrentPlotNo) ';']);
ORG.ExecuteLabTalk(['set %C -cl color(' colour ');']);
ORG.ExecuteLabTalk('set %C -w 0;');
ORG.ExecuteLabTalk('set %C -z 5;');
ORG.ExecuteLabTalk('set %C -k 2;');
ORG.ExecuteLabTalk(['set %C -c color(' colour ');']);
ORG.ExecuteLabTalk(['set %C -cf color(' colour ');']);

ORG.CreateGraphPage('PanPcat');

ORG.ExecuteLabTalk(['plotxy iy:=[' sheetName ']Sheet1!(1,12,13) plot:=202 ogl:=[' ORG.ActiveGraphName ']' num2str(ORG.CurrentLayer) '!;'] );
colour = 'green';
ORG.CurrentPlotNo = ORG.NoPlots + 1;
ORG.NoPlots = ORG.NoPlots + 1;
ORG.ActivatePage(ORG.ActiveGraphName);
ORG.ExecuteLabTalk(['layer.plot = ' num2str(ORG.CurrentPlotNo) ';']);
ORG.ExecuteLabTalk(['set %C -cl color(' colour ');']);
ORG.ExecuteLabTalk('set %C -w 0;');
ORG.ExecuteLabTalk('set %C -z 5;');
ORG.ExecuteLabTalk('set %C -k 2;');
ORG.ExecuteLabTalk(['set %C -c color(' colour ');']);
ORG.ExecuteLabTalk(['set %C -cf color(' colour ');']);

ORG.Disconnect;

end
% [Cap]
% [Res]
% [Cap]/[Res]
% GHK