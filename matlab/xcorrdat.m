%
%		xcorrdatfft.m
%		Douglas L. Jones
%		University of Illinois
%		June 13, 2009
%
%	xcorrdatfft.m: computes cross-correlation data over specified lag range
%			efficiently using FFTs
%
%	NOTES:
%		x must be at least length startoffset+corrlen-1
%		y must be at least length startoffset+corrlen+maxlag-1
%		startoffset must be > minlag
%
%
function [delay,maxcorr,xpower,ypower,xcorrout] = xcorrdat(x,y,startoffset,corrlen,minlag,maxlag);

%   SETUP

%   set up vectors
ylong = y(startoffset+minlag:startoffset+corrlen-1+maxlag);
ylen = length(ylong);
xpad = zeros(size(ylong));
xpad(1:corrlen) = x(startoffset:startoffset+corrlen-1);

corrxy = real(ifft(fft(ylong).*conj(fft(xpad))));	% compute cross-correlation using FFTs

%  compute powers in sliding window of y
ypowers = zeros(1,maxlag-minlag+1);
ypowers(1) = y(startoffset+minlag:startoffset+minlag+corrlen-1)*y(startoffset+minlag:startoffset+minlag+corrlen-1)';

for ii=2:maxlag-minlag+1,
   ypowers(ii) = ypowers(ii-1) - ylong(ii-1)*ylong(ii-1) + ylong(ii+corrlen-1)*ylong(ii+corrlen-1);
end

%  compute sliding normalized correlation coefficient
xcorrout = zeros(1,maxlag-minlag+1);

xpower = x(startoffset:startoffset+corrlen-1)*x(startoffset:startoffset+corrlen-1)';

xcorrout = corrxy(1:maxlag-minlag+1)./sqrt(xpower*ypowers + eps);

%jjj = 1;
%for iii=minlag:maxlag,
%  power2 = y(startoffset-iii:startoffset+corrlen-iii-1)*y(startoffset-iii:startoffset+corrlen-iii-1)';
%  xcorrout(jjj) = x(startoffset:startoffset+corrlen-1)*y(startoffset-iii:startoffset+corrlen-iii-1)'/(sqrt(xpower*ypower)+eps);
%  jjj = jjj + 1;
%end

[maxcorr,delay] = max(xcorrout);
delay = delay + minlag - 1;
xpower = x(startoffset:startoffset+corrlen-1)*x(startoffset:startoffset+corrlen-1)'/corrlen;
ypower = y(startoffset+delay:startoffset+delay+corrlen-1)*y(startoffset+delay:startoffset+delay+corrlen-1)'/corrlen;

%   DONE!

