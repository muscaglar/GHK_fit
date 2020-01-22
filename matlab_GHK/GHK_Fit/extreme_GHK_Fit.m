function [params,output,params_nlin] = extreme_GHK_Fit(cIDs,saltName,polarity,dirName,plotName,sheetName,Origin)

%eg: (hfcl4_10mM,'HfCl4','Anion');

R = 8.3144621; % J / mol K
F_c = 96485.3399; %C / mol
T = 300; %K
U = R*T/F_c;

switch polarity
    case {'cation'}
        beta = [100,1];
    case {'anion'}
        beta = [0.1,1000];
end

maxResConc = 10000000;
minResConc = 0.00000001;

if strcmp(saltName,'ZrCl4')
  %  maxResConc = 0.9;
end
% if strcmp(saltName,'KCl')
%     maxResConc = 1.1;
%     minResConc = 0.0005;
% end

options = optimset('FinDiffRelStep',1e-10,'Algorithm','trust-region-reflective','Display','on','MaxIter',1e4,'MaxFunEvals',1e4,'FunValCheck','on','TolFun',1e-100,'TolX',1e-100);

switch saltName
    case {'ALaCl3','CeCl3','LaCl3'}
        F = @(P,xdata)(U).*log((243.*xdata.^3.*(P(1)./P(2)).^2 + 27.*xdata.^3.*(P(1)./P(2)) + 2.*xdata.^3 + 729.*xdata.^2.*(P(1)./P(2)).^3 + 108.*xdata.^2.*(P(1)./P(2)) + 3.*xdata.^2 + (4.*(-9.*xdata.^2.*(P(1)./P(2)) - xdata.^2 + 9.*xdata.*(P(1)./P(2)) - xdata + 2).^3 + (243.*xdata.^3.*(P(1)./P(2)).^2 + 27.*xdata.^3.*(P(1)./P(2)) + 2.*xdata.^3 + 729.*xdata.^2.*(P(1)./P(2)).^3 + 108.*xdata.^2.*(P(1)./P(2)) + 3.*xdata.^2 + 486.*xdata.*(P(1)./P(2)).^2 + 27.*xdata.*(P(1)./P(2)) + 15.*xdata + 81.*(P(1)./P(2)) + 7).^2).^(1./2) + 486.*xdata.*(P(1)./P(2)).^2 + 27.*xdata.*(P(1)./P(2)) + 15.*xdata + 81.*(P(1)./P(2)) + 7).^(1./3)./(3.*2.^(1./3).*(3.*xdata.*(P(1)./P(2)) + 1)) - (2.^(1./3).*(-9.*xdata.^2.*(P(1)./P(2)) - xdata.^2 + 9.*xdata.*(P(1)./P(2)) - xdata + 2))./(3.*(3.*xdata.*(P(1)./P(2)) + 1).*(243.*xdata.^3.*(P(1)./P(2)).^2 + 27.*xdata.^3.*(P(1)./P(2)) + 2.*xdata.^3 + 729.*xdata.^2.*(P(1)./P(2)).^3 + 108.*xdata.^2.*(P(1)./P(2)) + 3.*xdata.^2 + sqrt(4.*(-9.*xdata.^2.*(P(1)./P(2)) - xdata.^2 + 9.*xdata.*(P(1)./P(2)) - xdata + 2).^3 + (243.*xdata.^3.*(P(1)./P(2)).^2 + 27.*xdata.^3.*(P(1)./P(2)) + 2.*xdata.^3 + 729.*xdata.^2.*(P(1)./P(2)).^3 + 108.*xdata.^2.*(P(1)./P(2)) + 3.*xdata.^2 + 486.*xdata.*(P(1)./P(2)).^2 + 27.*xdata.*(P(1)./P(2)) + 15.*xdata+ 81.*(P(1)./P(2))+ 7).^2) + 486.*xdata.*(P(1)./P(2)).^2 + 27.*xdata.*(P(1)./P(2)) + 15.*xdata + 81.*(P(1)./P(2))+ 7).^(1./3)) - (1 - xdata)./(3.*(3.*xdata.*(P(1)./P(2)) + 1)))+ (10E20*(P(1)./P(2) < 0)^2);
        if strcmp(polarity,'cation')
            prefix = -3;
        elseif strcmp(polarity,'anion')
            prefix = 1;
        end
    case   {'BKCl','KCl','LiCl'}
        F = @(P,xdata)(U).*log( ((P(1)./P(2)) + xdata) ./ ((P(1)./P(2)).*xdata + 1) )+ (10E20*(P(1)./P(2) < 0)^2);
        if strcmp(polarity,'cation')
            prefix = -1;
        elseif strcmp(polarity,'anion')
            prefix = 1;
        end
    case   {'MgCl2'}
        F = @(P,xdata)(U).*log((sqrt(8.*xdata.^2.*(P(1)./P(2)) + xdata.^2 + 16.*xdata.*(P(1)./P(2)).^2 + 2.*xdata + 8.*(P(1)./P(2)) + 1) + xdata - 1)./(2.*(2.*xdata.*(P(1)./P(2)) + 1)))+ (10E20*(P(1)./P(2) < 0)^2);
        if strcmp(polarity,'cation')
            prefix = -2;
        elseif strcmp(polarity,'anion')
            prefix = 1;
        end
    case   {'HfCl4','ZrCl4'}
        F = @(P,xdata)tetra_fit(xdata,P);
        if strcmp(polarity,'cation')
            prefix = -4;
        elseif strcmp(polarity,'anion')
            prefix = 1;
        end
    case   {'K3PO4'}
        F=@(P,xdata)(U).*log((7.*xdata.^3.*(P(1)./P(2)).^3+81.*xdata.^3.*(P(1)./P(2)).^2+15.*xdata.^2.*(P(1)./P(2)).^3+27.*xdata.^2.*(P(1)./P(2)).^2+486.*xdata.^2.*(P(1)./P(2))+((7.*xdata.^3.*(P(1)./P(2)).^3+81.*xdata.^3.*(P(1)./P(2)).^2+15.*xdata.^2.*(P(1)./P(2)).^3+27.*xdata.^2.*(P(1)./P(2)).^2+486.*xdata.^2.*(P(1)./P(2))+3.*xdata.*(P(1)./P(2)).^3+108.*xdata.*(P(1)./P(2)).^2+729.*xdata+2.*(P(1)./P(2)).^3+27.*(P(1)./P(2)).^2+243.*(P(1)./P(2))).^2+4.*(3.*(xdata-1).*(P(1)./P(2)).*(xdata.*(P(1)./P(2))+3)-(xdata-1).^2.*(P(1)./P(2)).^2).^3).^(1./2)+3.*xdata.*(P(1)./P(2)).^3+108.*xdata.*(P(1)./P(2)).^2+729.*xdata+2.*(P(1)./P(2)).^3+27.*(P(1)./P(2)).^2+243.*(P(1)./P(2))).^(1./3)./(3.*2.^(1./3).*(xdata.*(P(1)./P(2))+3))-(2.^(1./3).*(3.*(xdata-1).*(P(1)./P(2)).*(xdata.*(P(1)./P(2))+3)-(xdata-1).^2.*(P(1)./P(2)).^2))./(3.*(xdata.*(P(1)./P(2))+3).*(7.*xdata.^3.*(P(1)./P(2)).^3+81.*xdata.^3.*(P(1)./P(2)).^2+15.*xdata.^2.*(P(1)./P(2)).^3+27.*xdata.^2.*(P(1)./P(2)).^2+486.*xdata.^2.*(P(1)./P(2))+((7.*xdata.^3.*(P(1)./P(2)).^3+81.*xdata.^3.*(P(1)./P(2)).^2+15.*xdata.^2.*(P(1)./P(2)).^3+27.*xdata.^2.*(P(1)./P(2)).^2+486.*xdata.^2.*(P(1)./P(2))+3.*xdata.*(P(1)./P(2)).^3+108.*xdata.*(P(1)./P(2)).^2+729.*xdata+2.*(P(1)./P(2)).^3+27.*(P(1)./P(2)).^2+243.*(P(1)./P(2))).^2+4.*(3.*(xdata-1).*(P(1)./P(2)).*(xdata.*(P(1)./P(2))+3)-(xdata-1).^2.*(P(1)./P(2)).^2).^3).^(1./2)+3.*xdata.*(P(1)./P(2)).^3+108.*xdata.*(P(1)./P(2)).^2+729.*xdata+2.*(P(1)./P(2)).^3+27.*(P(1)./P(2)).^2+243.*(P(1)./P(2))).^(1./3))-((xdata-1).*(P(1)./P(2)))./(3.*(xdata.*(P(1)./P(2))+3)))+ (10E20*(P(1)./P(2) < 0)^2);
        if strcmp(polarity,'cation')
            prefix = -1;
        elseif strcmp(polarity,'anion')
            prefix = 3;
        end
    otherwise
        disp('Warning: Specify the salt or closest ideal model salt');
        return
end


if cIDs <= 0
    aResConcs = [1E-5, 1E-4, 1E-3, 1E-2, 1E-1, 1,10] ;
    aCapConcs = (zeros(1,size(aResConcs,2)) ) + cIDs ;
    x_conc = aResConcs ./ aCapConcs;
    y_v = prefix.*(R.*T ./ F_c) .* log(x_conc);
else
    [aResConcs,aCapConcs,aVoltageOffsets,~,errors, ~,~, ~] = ca_selectivity_alteredError(cIDs,0,0,minResConc,maxResConc);
    x_conc = aResConcs ./ aCapConcs;
    y_v = (-1.*aVoltageOffsets)*1E-3; %V
    std_mean_error_V = errors * 1E-3;
end

try
    [betaNew,RNew,JNew,CovBNew,MSENew,ErrorModelInfoNew] = nlinfit(x_conc,y_v,F,beta);
    [ci1, se1] = nlparci_edited(betaNew,RNew,'covar',CovBNew);
    [ci2, se2] = nlparci_edited(betaNew,RNew,'jacobian',JNew);
catch
    disp('nlinfit failed');
    betaNew = [0, 0];
    RNew =  [0, 0];
    JNew = [0, 0];
    CovBNew = [0, 0];
    MSENew = [0, 0];
    ErrorModelInfoNew = [0, 0];
    ci1 = [0, 0];
    ci2 = [0, 0];
end
[P,resnorm,resid,exitflag,output,lambda,J] = lsqcurvefit(F, beta,x_conc,y_v,[0,0],[],options);
[ci, se] = nlparci_edited(P,resid,'jacobian',J);

minBound = min(x_conc) / 10;
maxBound = max(x_conc) * 10;

x_model =[minBound:0.01:0.7,1.3:1:maxBound];

response = F(P,x_model);
response2 = F(betaNew,x_model);

output = catpad(2,x_conc',y_v',std_mean_error_V',x_model',response',response2','padval',0);

%P_1_error P_2_error

params = [aCapConcs(1),P(1),[0],P(2),[0],P(1)/P(2),P(2)/P(1)]; %,resnorm,resid,exitflag,output,lambda,J,ci];
params_nlin = [aCapConcs(1),betaNew(1),[0],betaNew(2),[0],betaNew(1)/betaNew(2),betaNew(2)/betaNew(1)];

if(Origin)
    extreme_GHK_Output_Origin(output, plotName,sheetName,dirName);
    extreme_GHK_Params_Origin(params,[sheetName 'P'],dirName);
    extreme_GHK_Params_Origin(params_nlin,[sheetName 'PA'],dirName);
end
    function yout = tetra_fit(xdata, B)
        yout = zeros(size(xdata));
        opt = optimset('display','off','MaxFunEvals', 100000);
        step_size = 0.001;
        for i=1:length(xdata)
            yout(i) = fsolve(@(V)(((exp((U.*V))+exp((2.*F_c.*V)./(R.*T))+exp((3.*F_c.*V)./(R.*T))+1).*(exp((F_c.*V)./(R.*T))-xdata(i)))./(4-4.*xdata(i).*exp((4.*F_c.*V)./(R.*T))))-(B(1)./B(2))+(10E20*(B(1)./B(2) < 0)^2),step_size,opt);
        end
    end

end

