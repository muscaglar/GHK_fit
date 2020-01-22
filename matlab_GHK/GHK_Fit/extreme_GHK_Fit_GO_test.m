function extreme_GHK_Fit_GO_test(saltName,polarity,varargin)


ORG = Matlab2OriginPlot();
ORG.cd_TopLevel();
dirName = ['GHK_response_' saltName] ;
ORG.mkdir(dirName);
ORG.cd(dirName);
% subDir = [saltName '_Anion'];
% ORG.mkdir(subDir);
ORG.Disconnect;

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
    [params,~] = extreme_GHK_Fit_test(varargin{iArg},saltName,polarity,dirName,plotName,sheetName,1);
    paramsTotal = [paramsTotal; params];
    
    close all;
    iArg = iArg + 2;
end
sheetName = [dirName 'params'];
extreme_GHK_Output_Origin_nlnfit(paramsTotal, sheetName,dirName);%,subDir);
end