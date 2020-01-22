function [fitresult, gof] = fit_di(x,y,doPlot,initial)

[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( '(0.0259).*log((sqrt(8.*x.^2.*(C./A) + x.^2 + 16.*x.*(C./A).^2 + 2.*x + 8.*(C./A) + 1) + x - 1)./(2.*(2.*x.*(C./A) + 1)))+ (10E20*(C./A < 0)^2);', 'independent', 'x', 'dependent', 'y' );
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


