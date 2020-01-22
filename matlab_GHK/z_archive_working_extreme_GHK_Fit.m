function [params,output] = z_archive_working_extreme_GHK_Fit(cIDs,ion,polarity)

R = 8.3144621; % J / mol K
F_c = 96485.3399; %C / mol
T = 300; %K
U = R*T/F_c;

beta = 10;

options = optimset('FinDiffRelStep',1e1,'Algorithm','trust-region-reflective','Display','on','MaxIter',1e100,'MaxFunEvals',1e100,'FunValCheck','on','TolFun',1e-100,'TolX',1e-100);

switch ion
    case {'ALaCl3','CeCl3','LaCl3'}
        if strcmp(polarity,'cation')
            prefix = -3;
            F = @(P,xdata)(U).*log((243.*xdata.^3.*P.^2 + 27.*xdata.^3.*P + 2.*xdata.^3 + 729.*xdata.^2.*P.^3 + 108.*xdata.^2.*P + 3.*xdata.^2 + (4.*(-9.*xdata.^2.*P - xdata.^2 + 9.*xdata.*P - xdata + 2).^3 + (243.*xdata.^3.*P.^2 + 27.*xdata.^3.*P + 2.*xdata.^3 + 729.*xdata.^2.*P.^3 + 108.*xdata.^2.*P + 3.*xdata.^2 + 486.*xdata.*P.^2 + 27.*xdata.*P + 15.*xdata + 81.*P + 7).^2).^(1./2) + 486.*xdata.*P.^2 + 27.*xdata.*P + 15.*xdata + 81.*P + 7).^(1./3)./(3.*2.^(1./3).*(3.*xdata.*P + 1)) - (2.^(1./3).*(-9.*xdata.^2.*P - xdata.^2 + 9.*xdata.*P - xdata + 2))./(3.*(3.*xdata.*P + 1).*(243.*xdata.^3.*P.^2 + 27.*xdata.^3.*P + 2.*xdata.^3 + 729.*xdata.^2.*P.^3 + 108.*xdata.^2.*P + 3.*xdata.^2 + sqrt(4.*(-9.*xdata.^2.*P - xdata.^2 + 9.*xdata.*P - xdata + 2).^3 + (243.*xdata.^3.*P.^2 + 27.*xdata.^3.*P + 2.*xdata.^3 + 729.*xdata.^2.*P.^3 + 108.*xdata.^2.*P + 3.*xdata.^2 + 486.*xdata.*P.^2 + 27.*xdata.*P + 15.*xdata+ 81.*P+ 7).^2) + 486.*xdata.*P.^2 + 27.*xdata.*P + 15.*xdata + 81.*P+ 7).^(1./3)) - (1 - xdata)./(3.*(3.*xdata.*P + 1)));
        elseif strcmp(polarity,'anion')
            prefix = 1;
            F = @(P,xdata)(U).*log((243.*xdata.^3.*(1/P).^2 + 27.*xdata.^3.*(1/P) + 2.*xdata.^3 + 729.*xdata.^2.*(1/P).^3 + 108.*xdata.^2.*(1/P) + 3.*xdata.^2 + (4.*(-9.*xdata.^2.*(1/P) - xdata.^2 + 9.*xdata.*(1/P) - xdata + 2).^3 + (243.*xdata.^3.*(1/P).^2 + 27.*xdata.^3.*(1/P) + 2.*xdata.^3 + 729.*xdata.^2.*(1/P).^3 + 108.*xdata.^2.*(1/P) + 3.*xdata.^2 + 486.*xdata.*(1/P).^2 + 27.*xdata.*(1/P) + 15.*xdata + 81.*(1/P) + 7).^2).^(1./2) + 486.*xdata.*(1/P).^2 + 27.*xdata.*(1/P) + 15.*xdata + 81.*(1/P) + 7).^(1./3)./(3.*2.^(1./3).*(3.*xdata.*(1/P) + 1)) - (2.^(1./3).*(-9.*xdata.^2.*(1/P) - xdata.^2 + 9.*xdata.*(1/P) - xdata + 2))./(3.*(3.*xdata.*(1/P) + 1).*(243.*xdata.^3.*(1/P).^2 + 27.*xdata.^3.*(1/P) + 2.*xdata.^3 + 729.*xdata.^2.*(1/P).^3 + 108.*xdata.^2.*(1/P) + 3.*xdata.^2 + sqrt(4.*(-9.*xdata.^2.*(1/P) - xdata.^2 + 9.*xdata.*(1/P) - xdata + 2).^3 + (243.*xdata.^3.*(1/P).^2 + 27.*xdata.^3.*(1/P) + 2.*xdata.^3 + 729.*xdata.^2.*(1/P).^3 + 108.*xdata.^2.*(1/P) + 3.*xdata.^2 + 486.*xdata.*(1/P).^2 + 27.*xdata.*(1/P) + 15.*xdata+ 81.*(1/P)+ 7).^2) + 486.*xdata.*(1/P).^2 + 27.*xdata.*(1/P) + 15.*xdata + 81.*(1/P)+ 7).^(1./3)) - (1 - xdata)./(3.*(3.*xdata.*(1/P) + 1)));
        end
    case   {'BKCl','KCl','LiCl'}
        if strcmp(polarity,'cation')
            prefix = -1;
            F = @(P,xdata)(U).*log( ((P(1)./P(2)) + xdata) ./ ((P(1)./P(2)).*xdata + 1) );
        elseif strcmp(polarity,'anion')
            prefix = 1;
            F = @(P,xdata)(U).*log( ((P(2)./P(1)) + xdata) ./ ((P(2)./P(1)).*xdata + 1) );
        end
    case   {'MgCl2'}
        if strcmp(polarity,'cation')
            prefix = -2;
            F = @(P,xdata)(U).*log((sqrt(8.*xdata.^2.*P + xdata.^2 + 16.*xdata.*P.^2 + 2.*xdata + 8.*P + 1) + xdata - 1)./(2.*(2.*xdata.*P + 1)));
        elseif strcmp(polarity,'anion')
            prefix = 1;
            F = @(P,xdata)(U).*log((sqrt(8.*xdata.^2.*(1/P) + xdata.^2 + 16.*xdata.*(1/P).^2 + 2.*xdata + 8.*(1/P) + 1) + xdata - 1)./(2.*(2.*xdata.*(1/P) + 1)));
        end
    case   {'HfCl4','ZrCl4'}
        if strcmp(polarity,'cation')
            prefix = -4;
            F = @(P,xdata)tetra_fit_pcat_an(xdata,P);
        elseif strcmp(polarity,'anion')
            prefix = 1;
            F = @(P,xdata)tetra_fit_pan_cat(xdata,P);
        end
    case   {'K3PO4'}
        if strcmp(polarity,'cation')
            prefix = -1;
            F=@(P,xdata)(U).*log((7.*xdata.^3.*P.^3+81.*xdata.^3.*P.^2+15.*xdata.^2.*P.^3+27.*xdata.^2.*P.^2+486.*xdata.^2.*P+((7.*xdata.^3.*P.^3+81.*xdata.^3.*P.^2+15.*xdata.^2.*P.^3+27.*xdata.^2.*P.^2+486.*xdata.^2.*P+3.*xdata.*P.^3+108.*xdata.*P.^2+729.*xdata+2.*P.^3+27.*P.^2+243.*P).^2+4.*(3.*(xdata-1).*P.*(xdata.*P+3)-(xdata-1).^2.*P.^2).^3).^(1./2)+3.*xdata.*P.^3+108.*xdata.*P.^2+729.*xdata+2.*P.^3+27.*P.^2+243.*P).^(1./3)./(3.*2.^(1./3).*(xdata.*P+3))-(2.^(1./3).*(3.*(xdata-1).*P.*(xdata.*P+3)-(xdata-1).^2.*P.^2))./(3.*(xdata.*P+3).*(7.*xdata.^3.*P.^3+81.*xdata.^3.*P.^2+15.*xdata.^2.*P.^3+27.*xdata.^2.*P.^2+486.*xdata.^2.*P+((7.*xdata.^3.*P.^3+81.*xdata.^3.*P.^2+15.*xdata.^2.*P.^3+27.*xdata.^2.*P.^2+486.*xdata.^2.*P+3.*xdata.*P.^3+108.*xdata.*P.^2+729.*xdata+2.*P.^3+27.*P.^2+243.*P).^2+4.*(3.*(xdata-1).*P.*(xdata.*P+3)-(xdata-1).^2.*P.^2).^3).^(1./2)+3.*xdata.*P.^3+108.*xdata.*P.^2+729.*xdata+2.*P.^3+27.*P.^2+243.*P).^(1./3))-((xdata-1).*P)./(3.*(xdata.*P+3)));
        elseif strcmp(polarity,'anion')
            prefix = 3;
            F=@(P,xdata)(U).*log((7.*xdata.^3.*(1/P).^3+81.*xdata.^3.*(1/P).^2+15.*xdata.^2.*(1/P).^3+27.*xdata.^2.*(1/P).^2+486.*xdata.^2.*(1/P)+((7.*xdata.^3.*(1/P).^3+81.*xdata.^3.*(1/P).^2+15.*xdata.^2.*(1/P).^3+27.*xdata.^2.*(1/P).^2+486.*xdata.^2.*(1/P)+3.*xdata.*(1/P).^3+108.*xdata.*(1/P).^2+729.*xdata+2.*(1/P).^3+27.*(1/P).^2+243.*(1/P)).^2+4.*(3.*(xdata-1).*(1/P).*(xdata.*(1/P)+3)-(xdata-1).^2.*(1/P).^2).^3).^(1./2)+3.*xdata.*(1/P).^3+108.*xdata.*(1/P).^2+729.*xdata+2.*(1/P).^3+27.*(1/P).^2+243.*(1/P)).^(1./3)./(3.*2.^(1./3).*(xdata.*(1/P)+3))-(2.^(1./3).*(3.*(xdata-1).*(1/P).*(xdata.*(1/P)+3)-(xdata-1).^2.*(1/P).^2))./(3.*(xdata.*(1/P)+3).*(7.*xdata.^3.*(1/P).^3+81.*xdata.^3.*(1/P).^2+15.*xdata.^2.*(1/P).^3+27.*xdata.^2.*(1/P).^2+486.*xdata.^2.*(1/P)+((7.*xdata.^3.*(1/P).^3+81.*xdata.^3.*(1/P).^2+15.*xdata.^2.*(1/P).^3+27.*xdata.^2.*(1/P).^2+486.*xdata.^2.*(1/P)+3.*xdata.*(1/P).^3+108.*xdata.*(1/P).^2+729.*xdata+2.*(1/P).^3+27.*(1/P).^2+243.*(1/P)).^2+4.*(3.*(xdata-1).*(1/P).*(xdata.*(1/P)+3)-(xdata-1).^2.*(1/P).^2).^3).^(1./2)+3.*xdata.*(1/P).^3+108.*xdata.*(1/P).^2+729.*xdata+2.*(1/P).^3+27.*(1/P).^2+243.*(1/P)).^(1./3))-((xdata-1).*(1/P))./(3.*(xdata.*(1/P)+3)));
        end
    otherwise
        disp('Warning: Specify the salt or closest ideal model salt');
        return
end


if cIDs <= 10
    aResConcs = [1E-5, 1E-4, 1E-3, 1E-2, 1E-1, 1,10] ;
    aCapConcs = (zeros(1,size(aResConcs,2)) ) + cIDs ;
    x_conc = aResConcs ./ aCapConcs;
    y_v = prefix.*(R.*T ./ F_c) .* log(x_conc);
else
    [aResConcs,aCapConcs,aVoltageOffsets,~,errors, ~,~, ~,std_mean_error_V] = ca_selectivity(cIDs,1,minResConc,maxResConc);
    x_conc = aResConcs ./ aCapConcs;
    y_v = (-1.*aVoltageOffsets)*1E-3; %V
    std_mean_error_V = std_mean_error_V * 1E-3;
end

[P,resnorm,resid,exitflag,output,lambda,J] = lsqcurvefit(F,[beta,beta],x_conc,y_v,[0, 0],[],options);
ci = nlparci(P,resid,'jacobian',J);
response = F(P,x_conc);
output = [x_conc',y_v',response'];
params = [P];%,resnorm,resid,exitflag,output,lambda,J,ci];

plot(x_conc,y_v,'ro');
hold on
plot(x_conc,response)
hold off
set(gca,'XScale','log')
figure;
plot(x_conc,y_v,'ro');
hold on
plot(x_conc,response-resid)
set(gca,'XScale','log')
hold off
figure;
plot(x_conc,y_v,'ro');
hold on
plot(x_conc,response+resid)
set(gca,'XScale','log')
hold off

function yout = tetra_fit_pcat_an(xdata, B)
    yout = zeros(size(xdata));
    opt = optimset('display','off','MaxFunEvals', 100000);
    step_size = 0.001;
    for i=1:length(xdata)
        yout(i) = fsolve(@(V)(((exp((U.*V))+exp((2.*F_c.*V)./(R.*T))+exp((3.*F_c.*V)./(R.*T))+1).*(exp((F_c.*V)./(R.*T))-xdata(i)))./(4-4.*xdata(i).*exp((4.*F_c.*V)./(R.*T))))-B,step_size,opt);
    end
end

function yout = tetra_fit_pan_cat(xdata, B)
    yout = zeros(size(xdata));
    opt = optimset('display','off','MaxFunEvals', 10000000000);
    step_size = 0.001;
    for i=1:length(xdata)
        yout(i) = fsolve(@(V)(((exp((F_c.*V)./(R.*T))+exp((2.*F_c.*V)./(R.*T))+exp((3.*F_c.*V)./(R.*T))+1).*(exp((F_c.*V)./(R.*T))-xdata(i)))./(4-4.*xdata(i).*exp((4.*F_c.*V)./(R.*T))))-(1./B),step_size,opt);
    % yout(i) = fsolve(@(V)(( ( ( exp(V.*U) + 1 ) .* (exp(2.*U.*V) + 1 ) .* (exp(U.*V) - xdata(i) ) ) ./ ( 4 - (4.*xdata(i).*exp(4.*U.*V) ) )  )  - (1./P)),0.001,opt);
    end
end

end

