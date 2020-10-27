function [Gs_sto_final] = modelreducerapp(Gs_s_m)

% Perform mode selection on LTI system
System = Gs_s_m; % Define System to reduce
UpperCutoffFrequency = 50;
LowerCutoffFrequency = 1;
 
% Create option set for freqsep command
Options = freqsepOptions();
 
% Select modes between lower and upper cutoff frequencies
Gs_sto_final = freqsep(System,UpperCutoffFrequency,Options);
[~,Gs_sto_final] = freqsep(Gs_sto_final,LowerCutoffFrequency,Options);
