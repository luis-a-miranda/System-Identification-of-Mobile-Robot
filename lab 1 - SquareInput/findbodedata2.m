function [gain, phase, w] = findbodedata2(chosen_ref, chosen_y, Ts, neglect)

% function [gain, phase, w] = findbodedata2(input, output, Ts, neglect)
% This function receives a response (output) to a sinusoid (input),
% both sampled with sampling time Ts (in s, default is 1 s),
% and calculates gain (dB), phase lag (º), and input and output frequency (w, in rad/s).
% The first neglect(1) points (default is 0) and the last neglect(2) points (default is 0) are not considered.
% Duarte Valério 2007-2013; see below for code incorporated from Matlab FileExchange

if nargin < 4 % by default no points are neglected anywhere
    neglect = [0 0];
end
if length(neglect) == 1 % by default no points are neglected at the end
    neglect = [neglect 0];
end
if length(chosen_ref) < sum(neglect)
    error('You cannot neglect all the data.')
end
if length(chosen_ref) ~= length(chosen_y)
    error('The two sinusoids must have the same lengths.')
end
if nargin < 3
    Ts = 1;
end

chosen_ref = chosen_ref(neglect(1)+1:end-neglect(2));
chosen_y = chosen_y(neglect(1)+1:end-neglect(2));
chosen_ref = chosen_ref - mean(chosen_ref);
chosen_y = chosen_y - mean(chosen_y);

gain = 20*log10( sqrt(mean(chosen_y.^2)) / sqrt(mean(chosen_ref.^2)) );

lengthSignal = length(chosen_ref); % number of samples
windowsize = 2^nextpow2(lengthSignal);
INPUT = fft(chosen_ref, windowsize) / lengthSignal; % Fourier Transform of input
[~, indexINPUT] = max( abs(INPUT(1 : windowsize/2+1)) );
phaseINPUT = angle(INPUT(indexINPUT)); % phase of the most significant frequency component of input
OUTPUT = fft(chosen_y, windowsize) / lengthSignal; % Fourier Transform of output
[~, indexOUTPUT] = max( abs(OUTPUT(1 : windowsize/2+1)) );
phaseOUTPUT = angle(OUTPUT(indexOUTPUT)); % phase of the most significant frequency component of output
phase = rad2deg(phaseOUTPUT - phaseINPUT);

f = (2*pi/Ts) * (0:windowsize/2) / windowsize;
w = [f(indexINPUT), f(indexOUTPUT)];

% ORIGINAL HEADER AND COPYRIGHT OF FILE GIVEPD.M FROM MATLAB FILE EXCHANGE, ADAPTED INTO LINES 34-41

% Filename: givepd.m
% Author: Shashank G. Sawant
% email: sgsawant@gmail.com
% Description: Given 2 sinusoidal signals of the 
% same frequency, the function gives the "phase difference" between the 
% 2 given signals 
% The phase difference is in RADIANS!!!!!!
% The output is limited to pd={-pi,pi}radians
% Time Stamp: 2010 October 19 2043hrs

% Copyright (c) 2010, Shashank Sawant
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.