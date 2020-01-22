function [] = extreme_GHK_Output_Origin_test(input,plotName,sheetName,folderDirectory)

ORG = Matlab2OriginPlot();
%output = catpad(2,x',y',x_model',fitresult(x_model),'padval',0);
ORG.cd_TopLevel();
ORG.cd(folderDirectory);

minIn = min(nonzeros(input(:,1)));
maxIn = max(nonzeros(input(:,1)));

minBound = minIn / 10;
maxBound = maxIn * 10;

sheetName = ORG.MatrixToOrigin(input, sheetName);

%1: Y
%2: D/C
%3: Y Error
%4: X

ORG.ExecuteLabTalk('wks.col1.lname$ = Concetration Ratio;');
ORG.ExecuteLabTalk('wks.col1.type = 4;');

ORG.ExecuteLabTalk('wks.col2.lname$ = Voltage Offset;');
ORG.ExecuteLabTalk('wks.col2.type = 1;');

ORG.ExecuteLabTalk('wks.col3.lname$ = Voltage Offset Error;');
ORG.ExecuteLabTalk('wks.col3.type = 3;');

ORG.ExecuteLabTalk('wks.col4.lname$ = Concetration Ratio;');
ORG.ExecuteLabTalk('wks.col4.type = 4;');

ORG.ExecuteLabTalk('wks.col5.lname$ = Model Response;');
ORG.ExecuteLabTalk('wks.col5.type = 1;');

ORG.ExecuteLabTalk(['win -a ' sheetName ' ']);

currentLayer = ORG.CurrentLayer;

ORG.CreateGraphPage(plotName);

ORG.ExecuteLabTalk(['plotxy iy:=[' sheetName ']Sheet1!(4,5) plot:=202 ogl:=[' ORG.ActiveGraphName ']' num2str(ORG.CurrentLayer) '!;'] );
colour = 'blue';
ORG.CurrentPlotNo = ORG.NoPlots + 1;
ORG.NoPlots = ORG.NoPlots + 1;
ORG.ActivatePage(ORG.ActiveGraphName);
ORG.ExecuteLabTalk(['layer.plot = ' num2str(ORG.CurrentPlotNo) ';']);
ORG.ExecuteLabTalk(['set %C -cl color(' colour ');']);
ORG.ExecuteLabTalk('set %C -w 1000;');
ORG.ExecuteLabTalk('set %C -z 0;');
ORG.ExecuteLabTalk(['set %C -c color(' colour ');']);
ORG.ExecuteLabTalk(['set %C -cf color(' colour ');']);

ORG.ExecuteLabTalk(['plotxy iy:=[' sheetName ']Sheet1!(1,2,3) plot:=202 ogl:=[' ORG.ActiveGraphName ']' num2str(ORG.CurrentLayer) '!;'] );
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

ORG.ExecuteLabTalk(['page.lname$= ' ORG.ActiveGraphName ';'] );
ORG.FormatGraph;
plotName = ORG.ActiveGraphName;
ORG.ActivatePage(ORG.ActiveGraphName)

ORG.ExecuteLabTalk('yl.text$ = \p125(Nernst Potential (mV))');
ORG.ExecuteLabTalk('xb.text$ = \p125(Concentration Ratio ([Res]/[Cap]))');

ORG.ExecuteLabTalk('layer.x.label.pt=24');
ORG.ExecuteLabTalk('layer.y.label.pt=24');

ORG.ExecuteLabTalk('layer.x.type = 2');

ORG.ExecuteLabTalk(['layer.x.from = ' num2str(minBound)]);
ORG.ExecuteLabTalk(['layer.x.to = ' num2str(maxBound)]);

ORG.ExecuteLabTalk('layer.x.label.numFormat = 2');

ORG.Disconnect;

end
% [Cap]
% [Res]
% [Cap]/[Res]
% GHK