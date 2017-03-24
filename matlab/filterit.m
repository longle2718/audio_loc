%
%		filterit.m
%		Douglas L. Jones
%		University of Illinois
%		March 23, 2007
%
%	filterit.m: simple FFT-based frequency-sampled filter
%
%   NOTES:
%
% - Fix it so it works right for minfreq = 0, maxfreq = fs/2
%
function [out] = filterit(x,minfreq,maxfreq,fs)

%   SETUP
transitionband = 100;	% taper filter down over 100 Hz

len = max(size(x));
channelct = min(size(x));

%   Initialize arrays

minfsamp = ceil(len*minfreq/fs)+1;
maxfsamp = floor(len*maxfreq/fs)+1;
minfsamptrans = ceil(len*(minfreq-transitionband)/fs) + 1;
if (minfsamptrans < 1), minfsamptrans = 1, end
maxfsamptrans = floor(len*(maxfreq+transitionband)/fs) + 1;
if (maxfsamptrans > len/2), maxfsamptrans = floor(len/2), end

filtvec = zeros(len,1);
filtvec(minfsamp:maxfsamp) = ones(maxfsamp-minfsamp+1,1);
filtvec(minfsamptrans:minfsamp) = 0.5*(cos(pi*[0:minfsamp-minfsamptrans]/(minfsamp-minfsamptrans))' - ones(minfsamp-minfsamptrans+1,1));
filtvec(maxfsamp:maxfsamptrans) = 0.5*(cos(pi*[0:maxfsamptrans-maxfsamp]/(maxfsamptrans-maxfsamp))' + ones(maxfsamptrans-maxfsamp+1,1));
%filtvec(minfsamp) = 0.5;
%filtvec(maxfsamp) = 0.5;

filtmat = 2*filtvec*ones(1,channelct);


%
%  APPLY ALGORITHM
%
out = real(ifft(fft(x,[],2).*filtmat',[],2));

%  DONE

