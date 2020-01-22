function [fitresult, gof] = fit_mono(x,y,doPlot,initial)

[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( '(0.0259).*log( ((C/A) + x) ./ ((C/A).*x + 1) )+ (10E20*(C/A < 0)^2)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.DiffMaxChange = 0.01;
opts.DiffMinChange = 1e-10;
opts.Display = 'Off';
opts.Lower = [0.0001 0.1];
opts.MaxFunEvals = 6000;
opts.MaxIter = 4000;
opts.StartPoint = initial;
opts.TolFun = 1e-10;
opts.TolX = 1e-10;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
%[myFit] = NonLinearModel.fit( xData, yData, ft,initial);
%[fitresult, gof] = fit( xData, yData, ft, opts );

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


