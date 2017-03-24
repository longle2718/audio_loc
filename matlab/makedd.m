%
%		makedd.m
%		Douglas L. Jones
%		University of Illinois
%		March 22, 2007
%
%	makedd.m: This code creates relative distances and delays for signals
%		to arbitrarily placed microphones
%
function [gains,delays] = makedd(miclocs,locgrid,fs);


%   SETUP

%   set parameters

csound = 343;       % m/sec

ans = size(miclocs);
channelct = ans(1);

%
%  COMPUTE IMPULSE RESPONSES TO EACH MIC FROM EACH LOCATION
%
ans = size(locgrid);
locct = ans(1);

gains = zeros(channelct,locct);
delays = zeros(channelct,locct);

for ii=1:locct,
    for iii=1:channelct,
        distance = sqrt((miclocs(iii,:) - locgrid(ii,:))*(miclocs(iii,:) - locgrid(ii,:))');
	delays(iii,ii) = distance*fs/csound;
        gains(iii,ii) = 1/distance;
    end
    gains(:,ii) = gains(:,ii)/sqrt(gains(:,ii)'*gains(:,ii));
end

delays = delays - ones(channelct,1)*mean(delays);

%   DONE!

