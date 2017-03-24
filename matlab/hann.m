%
%		hann.m
%		Douglas L. Jones
%		University of Illinois
%		October 25, 2008
%
%	hann.m: Who needs Matlab's signal processing toolbox?
%
function [w] = hann(len);


%   SETUP

%   set parameters

w = 0.5 - 0.5*cos(2*pi*[0:len-1]/(len-1))';

%   DONE!

