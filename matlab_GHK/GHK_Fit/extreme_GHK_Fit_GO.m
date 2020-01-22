function extreme_GHK_Fit_GO(saltName,varargin)


ORG = Matlab2OriginPlot();
ORG.cd_TopLevel();
dirName = ['GHK_response_' saltName] ;
ORG.mkdir(dirName);
ORG.cd(dirName);
% subDir = [saltName '_Anion'];
% ORG.mkdir(subDir);
ORG.Disconnect;

paramsTotal = [];
paramsAltTotal = [];

%extreme_GHK_Fit_GO('HfCl4',hfcl4_100mM,'100mM',hfcl4_10mM,'10mM',hfcl4_1mM,'1mM',hfcl4_1M,'1M')

iArg = 1;
while iArg < size(varargin,2)
    currentConc = varargin{iArg+1};
    sheetName = ['GHK_' saltName '_' currentConc];
    plotName = ['GHK_' saltName '_' currentConc];
    [params,~,paramsAlt] = extreme_GHK_Fit(varargin{iArg},saltName,'Anion',dirName,plotName,sheetName,1);
    paramsTotal = [paramsTotal; params];
    paramsAltTotal = [paramsAltTotal; paramsAlt];
    close all;
    iArg = iArg + 2;
end

sheetName = [dirName 'params'];
extreme_GHK_Params_Origin(paramsTotal, sheetName,dirName);%,subDir);

sheetName = [dirName 'paramsAlt'];
extreme_GHK_Params_Origin(paramsAltTotal, sheetName,dirName);%,subDir);
end