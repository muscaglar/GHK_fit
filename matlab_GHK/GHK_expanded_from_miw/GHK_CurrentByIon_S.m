function [ CurrentS ] = GHK_CurrentByIon_S( z, E , P, ConcI, Conc0 )
%GHK Summary of this function goes here
%   Detailed explanation goes here
PhysicalConsts

CurrentS = P .* z.^2 .* ( ( E .* F_Const.^2)  ./ (R_Const .* T_Const) ) .* ( ( ConcI - Conc0 .* exp(-1 .* GHK_Vs(z,E)) ) ./ ( 1 - exp(-1 .* GHK_Vs(z,E)) ) ); 

%CurrentS = ( ( ConcI - Conc0 .* exp(-1 .* GHK_Vs(z,E)) ) ./ ( 1 - exp(-1 .* GHK_Vs(z,E)) ) ); 

end

