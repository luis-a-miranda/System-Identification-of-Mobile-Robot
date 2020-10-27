function [Gs_30_final]=modelreducerapp30(Gs_m)

% Perform mode selection on LTI system
System = Gs_m; % Define System to reduce
UpperCutoffFrequency = 40;
LowerCutoffFrequency = 2.154434690031883;
 
% Create option set for freqsep command
Options = freqsepOptions();
 
% Select modes between lower and upper cutoff frequencies
Gs_30_final = freqsep(System,UpperCutoffFrequency,Options);
[~,Gs_30_final] = freqsep(Gs_30_final,LowerCutoffFrequency,Options);
 
