function [fitresult, gof] = fit_other_tri(x,y,doPlot,initial)

[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( '(0.0259).*log((7.*x.^3.*(C./A).^3+81.*x.^3.*(C./A).^2+15.*x.^2.*(C./A).^3+27.*x.^2.*(C./A).^2+486.*x.^2.*(C./A)+((7.*x.^3.*(C./A).^3+81.*x.^3.*(C./A).^2+15.*x.^2.*(C./A).^3+27.*x.^2.*(C./A).^2+486.*x.^2.*(C./A)+3.*x.*(C./A).^3+108.*x.*(C./A).^2+729.*x+2.*(C./A).^3+27.*(C./A).^2+243.*(C./A)).^2+4.*(3.*(x-1).*(C./A).*(x.*(C./A)+3)-(x-1).^2.*(C./A).^2).^3).^(1./2)+3.*x.*(C./A).^3+108.*x.*(C./A).^2+729.*x+2.*(C./A).^3+27.*(C./A).^2+243.*(C./A)).^(1./3)./(3.*2.^(1./3).*(x.*(C./A)+3))-(2.^(1./3).*(3.*(x-1).*(C./A).*(x.*(C./A)+3)-(x-1).^2.*(C./A).^2))./(3.*(x.*(C./A)+3).*(7.*x.^3.*(C./A).^3+81.*x.^3.*(C./A).^2+15.*x.^2.*(C./A).^3+27.*x.^2.*(C./A).^2+486.*x.^2.*(C./A)+((7.*x.^3.*(C./A).^3+81.*x.^3.*(C./A).^2+15.*x.^2.*(C./A).^3+27.*x.^2.*(C./A).^2+486.*x.^2.*(C./A)+3.*x.*(C./A).^3+108.*x.*(C./A).^2+729.*x+2.*(C./A).^3+27.*(C./A).^2+243.*(C./A)).^2+4.*(3.*(x-1).*(C./A).*(x.*(C./A)+3)-(x-1).^2.*(C./A).^2).^3).^(1./2)+3.*x.*(C./A).^3+108.*x.*(C./A).^2+729.*x+2.*(C./A).^3+27.*(C./A).^2+243.*(C./A)).^(1./3))-((x-1).*(C./A))./(3.*(x.*(C./A)+3)))+ (10E20*(C./A < 0)^2);', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0];
opts.StartPoint = initial;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

if(doPlot==1)
    % Plot fit with data.
    figure( 'Name', 'untitled fit 1' );
    plot( fitresult, xData, yData);
    % Label axes
    xlabel x
    ylabel y
    grid on
    
end
end


