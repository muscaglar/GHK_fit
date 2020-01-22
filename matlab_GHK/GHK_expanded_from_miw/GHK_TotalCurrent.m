
function [ I_Total, I_Components ] = GHK_TotalCurrent( z, E , P, ConcI, Conc0)
%GHK Total Current Input a vector of z and P  - for each then calc the
%current and sum

%Note this function expects and returns SI units. Ie Volts and Amps, not mV
%and nA

n = max(size(z));
I_Components = zeros(max(size(E)),n);
for i = 1:n
    I_Components(:,i) = GHK_CurrentByIon_S( z(i), E , P(i), ConcI, Conc0 );
end
I_Total = sum(I_Components,2);

DoPlot = 0;
if(DoPlot >0)
    ORG = Matlab2OriginPlot();
    ORG.HoldOff();
    ORG.PlotScatter(E,I_Components','Hi2', ORG.ColourPicker());
    ORG.ylabel('Current','nA');
    ORG.xlabel('Voltage','mV');
    ORG.HideActiveWkBk()
    %ORG.xComment('');
    ORG.HoldOn()
    ORG.PlotScatter(E,I_Total','Hi2', ORG.ColourPicker());
    
     hold off;
     plot(E,I_Components,' - ');
     hold all;
     plot(E,I_Total);

    xlabel('Voltage V');   
    ylabel('Current A')
    title(['GHK Plot']);
end

end

