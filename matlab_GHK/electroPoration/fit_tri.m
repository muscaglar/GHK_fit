function [fitresult, gof] = fit_tri(x,y,doPlot,initial)

[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( '(0.0259).*log((243.*x.^3.*(C/A).^2 + 27.*x.^3.*(C/A) + 2.*x.^3 + 729.*x.^2.*(C/A).^3 + 108.*x.^2.*(C/A) + 3.*x.^2 + (4.*(-9.*x.^2.*(C/A) - x.^2 + 9.*x.*(C/A) - x + 2).^3 + (243.*x.^3.*(C/A).^2 + 27.*x.^3.*(C/A) + 2.*x.^3 + 729.*x.^2.*(C/A).^3 + 108.*x.^2.*(C/A) + 3.*x.^2 + 486.*x.*(C/A).^2 + 27.*x.*(C/A) + 15.*x + 81.*(C/A) + 7).^2).^(1./2) + 486.*x.*(C/A).^2 + 27.*x.*(C/A) + 15.*x + 81.*(C/A) + 7).^(1./3)./(3.*2.^(1./3).*(3.*x.*(C/A) + 1)) - (2.^(1./3).*(-9.*x.^2.*(C/A) - x.^2 + 9.*x.*(C/A) - x + 2))./(3.*(3.*x.*(C/A) + 1).*(243.*x.^3.*(C/A).^2 + 27.*x.^3.*(C/A) + 2.*x.^3 + 729.*x.^2.*(C/A).^3 + 108.*x.^2.*(C/A) + 3.*x.^2 + sqrt(4.*(-9.*x.^2.*(C/A) - x.^2 + 9.*x.*(C/A) - x + 2).^3 + (243.*x.^3.*(C/A).^2 + 27.*x.^3.*(C/A) + 2.*x.^3 + 729.*x.^2.*(C/A).^3 + 108.*x.^2.*(C/A) + 3.*x.^2 + 486.*x.*(C/A).^2 + 27.*x.*(C/A) + 15.*x+ 81.*(C/A)+ 7).^2) + 486.*x.*(C/A).^2 + 27.*x.*(C/A) + 15.*x + 81.*(C/A)+ 7).^(1./3)) - (1 - x)./(3.*(3.*x.*(C/A) + 1)));', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares','Algorithm','Trust-Region','MaxFunEvals',1000,'MaxIter',1000);
opts.Display = 'Off';
opts.Lower = [1 1];
opts.StartPoint = initial;


% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

if(doPlot==1)
    % Plot fit with data.
    close all;
    figure( 'Name', 'untitled fit 1' );
    plot( fitresult, xData, yData);
    % Label axes
    xlabel x
    ylabel y
    grid on
    set(gca, 'XScale', 'log');
end
end%+ (10E20*(C/A < 0)^2


