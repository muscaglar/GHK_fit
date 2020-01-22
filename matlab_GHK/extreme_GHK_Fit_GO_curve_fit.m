function [paramsTotal] = extreme_GHK_Fit_GO_curve_fit(saltName,polarity,varargin)

doPlot = 1;

if(doPlot)
ORG = Matlab2OriginPlot();
ORG.cd_TopLevel();
dirName = ['GHK_response_' saltName] ;
ORG.mkdir(dirName);
ORG.cd(dirName);
ORG.Disconnect;
else
    dirName = 0;
end

paramsTotal = [];

if(strcmp(saltName,'KCl_HfCl4'))
    saltName = 'KCl';
end
if(strcmp(saltName,'KCl_LaCl3'))
    saltName = 'KCl';
end
iArg = 1;
while iArg < size(varargin,2)
    currentConc = varargin{iArg+1};
    sheetName = ['GHK_' saltName '_' currentConc];
    plotName = ['GHK_' saltName '_' currentConc];
    [params,~] = extreme_GHK_Fit_curve_fit(varargin{iArg},saltName,polarity,dirName,plotName,sheetName,doPlot);
    paramsTotal = [paramsTotal; params];
    
    close all;
    iArg = iArg + 2;
end
if(doPlot)
sheetName = [dirName 'params'];
extreme_GHK_Output_Origin_nlnfit(paramsTotal, sheetName,dirName);%,subDir);
end
end