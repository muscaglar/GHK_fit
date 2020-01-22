function [ Vs ] = GHK_Vs( z, E, F, T, R )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
PhysicalConsts
if nargin < 5
    %Gas constant.
    R = R_Const;
end
if nargin < 4
    %Set temperature as 300K
    T = T_Const;
end
if nargin < 3
    %Set Faraday http://en.wikipedia.org/wiki/Faraday_constant
    %Its the charge per Mol of electrodes
    F = F_Const ;
end
if nargin < 2
    %Set voltage to 1V if not defined
    E = 1;
end
if nargin < 1
    %Set z to 1 if not defined
    z = 1;
end
Vs = z .* ((E .* F) ./ ( R .* T ));

end

