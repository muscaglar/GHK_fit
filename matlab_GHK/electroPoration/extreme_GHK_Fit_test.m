function [params,output] = extreme_GHK_Fit_test(cIDs,saltName,polarity,dirName,plotName,sheetName,Origin)

%eg: (hfcl4_10mM,'HfCl4','Anion');

R = 8.3144621; % J / mol K
F_c = 96485.3399; %C / mol
T = 300; %K
U = R*T/F_c;
options = optimset('FinDiffRelStep',1e-10,'Algorithm','trust-region-reflective','Display','on','MaxIter',1e4,'MaxFunEvals',1e3,'FunValCheck','on','TolFun',1e-100,'TolX',1e-100);
flag = 0;

% P is Pcat/Pan
polarity = lower(polarity);

switch polarity
    case {'cation'}
        initial = [1000 1];
    case {'anion'}
        initial = [1 1000];
end

maxResConc = 10000000;
minResConc = 0.00000001;

switch saltName
    case {'ALaCl3','CeCl3','LaCl3'}
        if strcmp(polarity,'cation')
            prefix = -3;
        elseif strcmp(polarity,'anion')
            prefix = 1;
        end
    case   {'BKCl','KCl','LiCl'}
        if strcmp(polarity,'cation')
            prefix = -1;
        elseif strcmp(polarity,'anion')
            prefix = 1;
        end
    case   {'HfCl4','ZrCl4'}
        if strcmp(polarity,'cation')
            prefix = -4;
        elseif strcmp(polarity,'anion')
            prefix = 1;
        end
    case   {'MgCl2'}
        if strcmp(polarity,'cation')
            prefix = -2;
        elseif strcmp(polarity,'anion')
            prefix = 1;
        end
    case   {'K3PO4'}
        if strcmp(polarity,'cation')
            prefix = -1;
        elseif strcmp(polarity,'anion')
            prefix = 3;
        end
end

if cIDs <= 10
    aResConcs = [1E-5, 1E-4, 1E-3, 1E-2, 1E-1, 1,10] ;
    aCapConcs = (zeros(1,size(aResConcs,2)) ) + cIDs ;
    x = aResConcs ./ aCapConcs;
    y = prefix.*(R.*T ./ F_c) .* log(x);
else
    [aResConcs,aCapConcs,aVoltageOffsets,~,errors, ~,~, ~,std_mean_error_V] = ca_selectivity(cIDs,1,minResConc,maxResConc);
    x = aResConcs ./ aCapConcs;
    y = (-1.*aVoltageOffsets)*1E-3; %V
    std_mean_error_V = std_mean_error_V * 1E-3;
end

minBound = min(x) / 10;
maxBound = max(x) * 10;

x_model =[minBound:0.01:0.7,1.3:1:maxBound];

MODEL_aResConcs = [1E-5, 1E-4, 1E-3, 1E-2, 1E-1, 1,10] ;
MODEL_aCapConcs = (zeros(1,size(MODEL_aResConcs,2)) ) + aCapConcs(1) ;
MODEL_x = MODEL_aResConcs ./ MODEL_aCapConcs;
MODEL_y = prefix.*(R.*T ./ F_c) .* log(MODEL_x);

%%Curve Fit Method
%[fitresult, gof] = fit_mono(MODEL_x,MODEL_y,plot,initial);
%initial = [fitresult.C,fitresult.A];
%[fitresult, gof] = fit_mono(x,y,plot,initial);

switch saltName
    case {'ALaCl3','CeCl3','LaCl3'}
        F = @(P,xdata)(U).*log((243.*xdata.^3.*(P(1)./P(2)).^2 + 27.*xdata.^3.*(P(1)./P(2)) + 2.*xdata.^3 + 729.*xdata.^2.*(P(1)./P(2)).^3 + 108.*xdata.^2.*(P(1)./P(2)) + 3.*xdata.^2 + (4.*(-9.*xdata.^2.*(P(1)./P(2)) - xdata.^2 + 9.*xdata.*(P(1)./P(2)) - xdata + 2).^3 + (243.*xdata.^3.*(P(1)./P(2)).^2 + 27.*xdata.^3.*(P(1)./P(2)) + 2.*xdata.^3 + 729.*xdata.^2.*(P(1)./P(2)).^3 + 108.*xdata.^2.*(P(1)./P(2)) + 3.*xdata.^2 + 486.*xdata.*(P(1)./P(2)).^2 + 27.*xdata.*(P(1)./P(2)) + 15.*xdata + 81.*(P(1)./P(2)) + 7).^2).^(1./2) + 486.*xdata.*(P(1)./P(2)).^2 + 27.*xdata.*(P(1)./P(2)) + 15.*xdata + 81.*(P(1)./P(2)) + 7).^(1./3)./(3.*2.^(1./3).*(3.*xdata.*(P(1)./P(2)) + 1)) - (2.^(1./3).*(-9.*xdata.^2.*(P(1)./P(2)) - xdata.^2 + 9.*xdata.*(P(1)./P(2)) - xdata + 2))./(3.*(3.*xdata.*(P(1)./P(2)) + 1).*(243.*xdata.^3.*(P(1)./P(2)).^2 + 27.*xdata.^3.*(P(1)./P(2)) + 2.*xdata.^3 + 729.*xdata.^2.*(P(1)./P(2)).^3 + 108.*xdata.^2.*(P(1)./P(2)) + 3.*xdata.^2 + sqrt(4.*(-9.*xdata.^2.*(P(1)./P(2)) - xdata.^2 + 9.*xdata.*(P(1)./P(2)) - xdata + 2).^3 + (243.*xdata.^3.*(P(1)./P(2)).^2 + 27.*xdata.^3.*(P(1)./P(2)) + 2.*xdata.^3 + 729.*xdata.^2.*(P(1)./P(2)).^3 + 108.*xdata.^2.*(P(1)./P(2)) + 3.*xdata.^2 + 486.*xdata.*(P(1)./P(2)).^2 + 27.*xdata.*(P(1)./P(2)) + 15.*xdata+ 81.*(P(1)./P(2))+ 7).^2) + 486.*xdata.*(P(1)./P(2)).^2 + 27.*xdata.*(P(1)./P(2)) + 15.*xdata + 81.*(P(1)./P(2))+ 7).^(1./3)) - (1 - xdata)./(3.*(3.*xdata.*(P(1)./P(2)) + 1)))+ (10E20*(P(1)./P(2) < 0)^2);
        [mdl,error,pValue,tStat,P,model_P] = execute_Fit(MODEL_x,MODEL_y,x,y, F,initial,saltName);
        flag = 1;
    case   {'BKCl','KCl','LiCl'}
        F = @(P,xdata)(U).*log( ((P(1)./P(2)) + xdata) ./ ((P(1)./P(2)).*xdata + 1) )+ (10E20*(P(1)./P(2) < 0)^2);
        % Remove Zero and first data point
        %x(length(x):length(x)-1)=[];
        %y(length(y):length(y)-1)=[];
        [mdl,error,pValue,tStat,P,model_P] = execute_Fit(MODEL_x,MODEL_y,x,y, F,initial,saltName);
        flag = 1;
    case   {'HfCl4','ZrCl4'}
        F = @(P,xdata)tetra_fit(xdata,P);
        % Remove Zero and first data point
        %         zero_ind = isalmost(MODEL_y,0,1E-5);
        %         zero_ind(1) = 1;
        %         MODEL_x(zero_ind) = [];
        %         MODEL_y(zero_ind) = [];
        [mdl,error,pValue,tStat,P,model_P] = execute_Fit(MODEL_x,MODEL_y,x,y, F,initial,saltName);
        flag = 1;
    case   {'MgCl2'}
        F = @(P,xdata)(U).*log((sqrt(8.*xdata.^2.*(P(1)./P(2)) + xdata.^2 + 16.*xdata.*(P(1)./P(2)).^2 + 2.*xdata + 8.*(P(1)./P(2)) + 1) + xdata - 1)./(2.*(2.*xdata.*(P(1)./P(2)) + 1)))+ (10E20*(P(1)./P(2) < 0)^2);
        [mdl,error,pValue,tStat,P,model_P] = execute_Fit(MODEL_x,MODEL_y,x,y, F,initial,saltName);
        flag = 1;
    case   {'K3PO4'}
        F=@(P,xdata)(U).*log((7.*xdata.^3.*(P(1)./P(2)).^3+81.*xdata.^3.*(P(1)./P(2)).^2+15.*xdata.^2.*(P(1)./P(2)).^3+27.*xdata.^2.*(P(1)./P(2)).^2+486.*xdata.^2.*(P(1)./P(2))+((7.*xdata.^3.*(P(1)./P(2)).^3+81.*xdata.^3.*(P(1)./P(2)).^2+15.*xdata.^2.*(P(1)./P(2)).^3+27.*xdata.^2.*(P(1)./P(2)).^2+486.*xdata.^2.*(P(1)./P(2))+3.*xdata.*(P(1)./P(2)).^3+108.*xdata.*(P(1)./P(2)).^2+729.*xdata+2.*(P(1)./P(2)).^3+27.*(P(1)./P(2)).^2+243.*(P(1)./P(2))).^2+4.*(3.*(xdata-1).*(P(1)./P(2)).*(xdata.*(P(1)./P(2))+3)-(xdata-1).^2.*(P(1)./P(2)).^2).^3).^(1./2)+3.*xdata.*(P(1)./P(2)).^3+108.*xdata.*(P(1)./P(2)).^2+729.*xdata+2.*(P(1)./P(2)).^3+27.*(P(1)./P(2)).^2+243.*(P(1)./P(2))).^(1./3)./(3.*2.^(1./3).*(xdata.*(P(1)./P(2))+3))-(2.^(1./3).*(3.*(xdata-1).*(P(1)./P(2)).*(xdata.*(P(1)./P(2))+3)-(xdata-1).^2.*(P(1)./P(2)).^2))./(3.*(xdata.*(P(1)./P(2))+3).*(7.*xdata.^3.*(P(1)./P(2)).^3+81.*xdata.^3.*(P(1)./P(2)).^2+15.*xdata.^2.*(P(1)./P(2)).^3+27.*xdata.^2.*(P(1)./P(2)).^2+486.*xdata.^2.*(P(1)./P(2))+((7.*xdata.^3.*(P(1)./P(2)).^3+81.*xdata.^3.*(P(1)./P(2)).^2+15.*xdata.^2.*(P(1)./P(2)).^3+27.*xdata.^2.*(P(1)./P(2)).^2+486.*xdata.^2.*(P(1)./P(2))+3.*xdata.*(P(1)./P(2)).^3+108.*xdata.*(P(1)./P(2)).^2+729.*xdata+2.*(P(1)./P(2)).^3+27.*(P(1)./P(2)).^2+243.*(P(1)./P(2))).^2+4.*(3.*(xdata-1).*(P(1)./P(2)).*(xdata.*(P(1)./P(2))+3)-(xdata-1).^2.*(P(1)./P(2)).^2).^3).^(1./2)+3.*xdata.*(P(1)./P(2)).^3+108.*xdata.*(P(1)./P(2)).^2+729.*xdata+2.*(P(1)./P(2)).^3+27.*(P(1)./P(2)).^2+243.*(P(1)./P(2))).^(1./3))-((xdata-1).*(P(1)./P(2)))./(3.*(xdata.*(P(1)./P(2))+3)))+ (10E20*(P(1)./P(2) < 0)^2);
        [mdl,error,pValue,tStat,P,model_P] = execute_Fit(MODEL_x,MODEL_y,x,y, F,initial,saltName);
        flag = 1;
    otherwise
        disp('Warning: Specify the salt or closest ideal model salt');
        return
end


if(flag)
    output = catpad(2,x',y',std_mean_error_V',x_model',F(P',x_model)','padval',0);
    extreme_GHK_Output_Origin_test(output, plotName,sheetName,dirName);
    
    Pcat = P(1);
    PcatErr = error(1);
    Pan = P(2);
    PanErr = error(2);
    PcatPan = P(1)/P(2);
    PanPcat = P(2)/P(1);
    
    PcatanErr = (sqrt((PcatErr/Pcat)^2 + (PanErr/Pan)^2))*PcatPan;
    PancatErr = (sqrt((PcatErr/Pcat)^2 + (PanErr/Pan)^2))*PanPcat;
    
    params = catpad(2,aCapConcs(1),Pcat,PcatErr,pValue(1),tStat(1),Pan,PanErr,pValue(2),tStat(2),PcatPan,PcatanErr,PanPcat,PancatErr,'padval',0);
    
    %extreme_GHK_Output_Origin_nlnfit(params,[sheetName 'P'],dirName);
    flag = 0;
else
    output = catpad(2,x',y',std_mean_error_V',x_model',fitresult(x_model),'padval',0);
    extreme_GHK_Output_Origin_test(output, plotName,sheetName,dirName);
    params = [aCapConcs(1),fitresult.C,fitresult.A,gof.sse,gof.adjrsquare,gof.rmse]; %,resnorm,resid,exitflag,output,lambda,J,ci];
    extreme_GHK_Params_Origin_test(params,[sheetName 'P'],dirName);
end

    function yout = tetra_fit(xdata, B)
        
        R = 8.3144621; % J / mol K
        F_c = 96485.3399; %C / mol
        T = 300; %K
        U = R*T/F_c;
        yout = zeros(size(xdata));
        opt = optimset('display','off','MaxFunEvals', 100000);
        step_size = 0.001;
        for i=1:length(xdata)
            yout(i) = fsolve(@(V)(((exp((U.*V))+exp((2.*F_c.*V)./(R.*T))+exp((3.*F_c.*V)./(R.*T))+1).*(exp((F_c.*V)./(R.*T))-xdata(i)))./(4-4.*xdata(i).*exp((4.*F_c.*V)./(R.*T))))-(B(1)./B(2))+(10E200*(B(1)./B(2) < 0)^2),step_size,opt);
        end
    end

    function [mdl,error,pValue,tStat,P,mod_P] = execute_Fit(mod_x,mod_y,real_x,real_y, equation,Pstart,salt)
        flag_tetra = 0;
        input_x = 0;
        
        switch salt
            case {'ALaCl3','CeCl3','LaCl3'}
                func = @(x)(U).*log((243.*input_x.^3.*(x(1)./x(2)).^2 + 27.*input_x.^3.*(x(1)./x(2)) + 2.*input_x.^3 + 729.*input_x.^2.*(x(1)./x(2)).^3 + 108.*input_x.^2.*(x(1)./x(2)) + 3.*input_x.^2 + (4.*(-9.*input_x.^2.*(x(1)./x(2)) - input_x.^2 + 9.*input_x.*(x(1)./x(2)) - input_x + 2).^3 + (243.*input_x.^3.*(x(1)./x(2)).^2 + 27.*input_x.^3.*(x(1)./x(2)) + 2.*input_x.^3 + 729.*input_x.^2.*(x(1)./x(2)).^3 + 108.*input_x.^2.*(x(1)./x(2)) + 3.*input_x.^2 + 486.*input_x.*(x(1)./x(2)).^2 + 27.*input_x.*(x(1)./x(2)) + 15.*input_x + 81.*(x(1)./x(2)) + 7).^2).^(1./2) + 486.*input_x.*(x(1)./x(2)).^2 + 27.*input_x.*(x(1)./x(2)) + 15.*input_x + 81.*(x(1)./x(2)) + 7).^(1./3)./(3.*2.^(1./3).*(3.*input_x.*(x(1)./x(2)) + 1)) - (2.^(1./3).*(-9.*input_x.^2.*(x(1)./x(2)) - input_x.^2 + 9.*input_x.*(x(1)./x(2)) - input_x + 2))./(3.*(3.*input_x.*(x(1)./x(2)) + 1).*(243.*input_x.^3.*(x(1)./x(2)).^2 + 27.*input_x.^3.*(x(1)./x(2)) + 2.*input_x.^3 + 729.*input_x.^2.*(x(1)./x(2)).^3 + 108.*input_x.^2.*(x(1)./x(2)) + 3.*input_x.^2 + sqrt(4.*(-9.*input_x.^2.*(x(1)./x(2)) - input_x.^2 + 9.*input_x.*(x(1)./x(2)) - input_x + 2).^3 + (243.*input_x.^3.*(x(1)./x(2)).^2 + 27.*input_x.^3.*(x(1)./x(2)) + 2.*input_x.^3 + 729.*input_x.^2.*(x(1)./x(2)).^3 + 108.*input_x.^2.*(x(1)./x(2)) + 3.*input_x.^2 + 486.*input_x.*(x(1)./x(2)).^2 + 27.*input_x.*(x(1)./x(2)) + 15.*input_x+ 81.*(x(1)./x(2))+ 7).^2) + 486.*input_x.*(x(1)./x(2)).^2 + 27.*input_x.*(x(1)./x(2)) + 15.*input_x + 81.*(x(1)./x(2))+ 7).^(1./3)) - (1 - input_x)./(3.*(3.*input_x.*(x(1)./x(2)) + 1)))+ (10E20*(x(1)./x(2) < 0)^2);
            case   {'BKCl','KCl','LiCl'}
                func = @(x)(U).*log( ((x(1)./x(2)) + input_x) ./ ((x(1)./x(2)).*input_x + 1) )+ (10E20*(x(1)./x(2) < 0)^2);
            case   {'HfCl4','ZrCl4'}
                flag_tetra = 1;
            case   {'MgCl2'}
                func = @(x)(U).*log((sqrt(8.*input_x.^2.*(x(1)./x(2)) + input_x.^2 + 16.*input_x.*(x(1)./x(2)).^2 + 2.*input_x + 8.*(x(1)./x(2)) + 1) + input_x - 1)./(2.*(2.*input_x.*(x(1)./x(2)) + 1)))+ (10E20*(x(1)./x(2) < 0)^2);
            case   {'K3PO4'}
                func =@(x)(U).*log((7.*input_x.^3.*(x(1)./x(2)).^3+81.*input_x.^3.*(x(1)./x(2)).^2+15.*input_x.^2.*(x(1)./x(2)).^3+27.*input_x.^2.*(x(1)./x(2)).^2+486.*input_x.^2.*(x(1)./x(2))+((7.*input_x.^3.*(x(1)./x(2)).^3+81.*input_x.^3.*(x(1)./x(2)).^2+15.*input_x.^2.*(x(1)./x(2)).^3+27.*input_x.^2.*(x(1)./x(2)).^2+486.*input_x.^2.*(x(1)./x(2))+3.*input_x.*(x(1)./x(2)).^3+108.*input_x.*(x(1)./x(2)).^2+729.*input_x+2.*(x(1)./x(2)).^3+27.*(x(1)./x(2)).^2+243.*(x(1)./x(2))).^2+4.*(3.*(input_x-1).*(x(1)./x(2)).*(input_x.*(x(1)./x(2))+3)-(input_x-1).^2.*(x(1)./x(2)).^2).^3).^(1./2)+3.*input_x.*(x(1)./x(2)).^3+108.*input_x.*(x(1)./x(2)).^2+729.*input_x+2.*(x(1)./x(2)).^3+27.*(x(1)./x(2)).^2+243.*(x(1)./x(2))).^(1./3)./(3.*2.^(1./3).*(input_x.*(x(1)./x(2))+3))-(2.^(1./3).*(3.*(input_x-1).*(x(1)./x(2)).*(input_x.*(x(1)./x(2))+3)-(input_x-1).^2.*(x(1)./x(2)).^2))./(3.*(input_x.*(x(1)./x(2))+3).*(7.*input_x.^3.*(x(1)./x(2)).^3+81.*input_x.^3.*(x(1)./x(2)).^2+15.*input_x.^2.*(x(1)./x(2)).^3+27.*input_x.^2.*(x(1)./x(2)).^2+486.*input_x.^2.*(x(1)./x(2))+((7.*input_x.^3.*(x(1)./x(2)).^3+81.*input_x.^3.*(x(1)./x(2)).^2+15.*input_x.^2.*(x(1)./x(2)).^3+27.*input_x.^2.*(x(1)./x(2)).^2+486.*input_x.^2.*(x(1)./x(2))+3.*input_x.*(x(1)./x(2)).^3+108.*input_x.*(x(1)./x(2)).^2+729.*input_x+2.*(x(1)./x(2)).^3+27.*(x(1)./x(2)).^2+243.*(x(1)./x(2))).^2+4.*(3.*(input_x-1).*(x(1)./x(2)).*(input_x.*(x(1)./x(2))+3)-(input_x-1).^2.*(x(1)./x(2)).^2).^3).^(1./2)+3.*input_x.*(x(1)./x(2)).^3+108.*input_x.*(x(1)./x(2)).^2+729.*input_x+2.*(x(1)./x(2)).^3+27.*(x(1)./x(2)).^2+243.*(x(1)./x(2))).^(1./3))-((input_x-1).*(x(1)./x(2)))./(3.*(input_x.*(x(1)./x(2))+3)))+ (10E20*(x(1)./x(2) < 0)^2);
        end
        
        answers = [];
        opts = statset('Display','iter','TolFun',1e-10,'TolX',1e-10,'MaxIter',400);
        
        %'Perfect' Fit
        perfect_mdl = fitnlm(mod_x,mod_y,equation,Pstart,'Options',opts);
        mod_P = perfect_mdl.Coefficients.Estimate;
        if (flag_tetra)
            x1 = mod_P(1);
            x2 = mod_P(2);
            flag_tetra = 0;
        else
            for input_x = mod_x
                nvars = 2;    % Number of variables
                LB = [1 1];   % Lower bound
                [P_ga,fval] = ga(func,nvars,[],[],[],[],LB,[],[]);
                answers = [answers; P_ga];
            end
            
            x1 = mean(answers(:,1));
            x2 = mean(answers(:,2));
        end
        
        noiseSignal = randi([95 105],1,2)./100;
        
        %Real Fit
        try
            mdl = fitnlm(real_x,real_y,equation,[x1,x2].*noiseSignal,'Options',opts);
            P = abs(mdl.Coefficients.Estimate);
            error = mdl.Coefficients.SE;
            pValue = mdl.Coefficients.pValue;
            tStat = mdl.Coefficients.tStat;
        catch
            disp('crashed')
            salt
            P = [1;1];
            error = [1;1];
            pValue =[1;1];
            tStat = [1;1];
            mdl = [];
        end
        
        
%         if(P(1)<P(2))
%             P(2) = P(2)/P(1);
%             P(1) = P(1)/P(1);
%         else
%             P(1) = P(1)/P(2);
%             P(2) = P(2)/P(2);
%         end
        
        
        % close all;
        % figure;
        % plot(mod_x,mod_y);
        % set(gca, 'XScale', 'log');
        % hold on;
        % scatter(real_x,real_y);
        % plot(mod_x,equation(P,mod_x));
        
    end
end

