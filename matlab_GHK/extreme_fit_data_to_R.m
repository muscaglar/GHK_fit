function [filenames] = extreme_fit_data_to_R(saltName,varargin)
flag = 0;
if(strcmp(saltName,'KCl_HfCl4'))
    saltName = 'KCl';
    flag = 1;
end
if(strcmp(saltName,'KCl_LaCl3'))
    saltName = 'KCl';
    flag = 2;
end
iArg = 1;
index = 1;

while iArg < size(varargin,2)
    currentConc = varargin{iArg+1};
    if(flag == 1)
         filename = ['C:\Users\mc934\Dropbox (Team Caglar)\data\GHK_' saltName 'HfCl4_' currentConc '.csv'];
    elseif (flag == 2 )
         filename = ['C:\Users\mc934\Dropbox (Team Caglar)\data\GHK_' saltName 'LaCl3_' currentConc '.csv'];
    else
    filename = ['C:\Users\mc934\Dropbox (Team Caglar)\data\GHK_' saltName '_' currentConc '.csv'];
    end
    extreme_fit_data_to_csv(varargin{iArg},0.000001,10000,filename);
    close all;
    filenames(index) = string(filename);
    iArg = iArg + 2;
    index =  index + 1;
end

    function [] = extreme_fit_data_to_csv(cIDs,minResConc,maxResConc,fName)
        [aResConcs,aCapConcs,aVoltageOffsets,~,errors,  ~, ~,  ~] = ca_selectivity_alteredError(cIDs,0,0,minResConc,maxResConc);
        std_mean_error_V = errors;
        x = aResConcs ./ aCapConcs;
        y = (-1.*aVoltageOffsets)*1E-3; %V
        error = std_mean_error_V * 1E-3;
        
        
        textHeader = "x,y,error";
        fid = fopen(fName,'w');
        fprintf(fid,'%s\n',textHeader)
        fclose(fid)
        dlmwrite(fName,[x',y',error'],'-append');
        
    end

end