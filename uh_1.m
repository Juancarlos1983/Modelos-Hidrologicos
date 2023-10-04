function [ u ] = uh_1( delay )
%uh_1_half Unit Hydrograph [days] with half a bell curve. GR4J-based
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
%   Inputs
%   d_base  - time base of routing delay [d]
%
%   Unit hydrograph spreads the input volume over a time period x4.
%   Percentage of input returned only increases. 
%   I.e. d_base = 3.8 [days], delta_t = 1: 
%   UH(1) = 0.04  [% of inflow]
%   UH(2) = 0.17
%   UH(3) = 0.35
%   UH(4) = 0.45

%%TIME STEP SIZE
if delay == 0
    delay = 1;                      % any value below t = 1 means no delay, but zero leads to problems
end                                           
tt = 1:ceil(delay);                 % Time series for which we need UH ordinates [days] 

%%EMPTIES
SH = zeros(1,length(tt) + 1); SH(1) = 0;
UH = zeros(1,length(tt));

%%UNIT HYDROGRAPH
for t = tt
    if t < delay
        SH(t + 1) = (t./delay).^(5./2);
    elseif t >= delay
        SH(t + 1) = 1;
    end    
    UH(t) = SH(t + 1) - SH(t);
end

%%DISPERSE VOLUME
u = UH;
end