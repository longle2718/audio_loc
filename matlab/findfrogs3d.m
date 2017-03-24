%
%		findfrogs3D.m
%		Douglas L. Jones
%		University of Illinois
%		September 6, 2008
%
%	findfrogs3D.m: Uses cross-correlation delays in frequency bands
%	  to locate single acoustic sources, shows where they are
%	  and when they're active, estimates steering vectors, and recovers them.
%    Extends findfrogs.m for frogs out of plane
%    - Extended on July 19, 2009 to compute all pairwise xcorrs and delays
%
%	USER NOTES
%	- Figure 10 plots the entire raw data record from startsamp to endsamp
%	  "cyan" corresponds to channel 1, "blue" to channel 2, "green" to channel 3, and "red" to channel 4
%	-Figure 11 - 10 + number_of_bands plots the bandpass data in each band
%	- Figure 1 maps the estimated source location on an x-y grid defined by
%	  parameters xmin,xmax,xinc,ymin,ymax,yinc, plotting 1/delay_error for
%	  each position in the grid, where delay_error is the difference between
%	  the measured relative delays between channels and the relative delays
%	  that should ideally be observed at that location.
%	- Figure 2 maps delay_error at each location (SMALL is good here),
%	  but gives an all-dark-blue map if the signal energy is considered
%	  too low to be real signal (this level controlled by parameter powersratiothresh),
%	  and an all-green map if the algorithm thinks the data is noisy or
%	  unreliable due to multiple sources in the analysis interval
%	  (as determined by the correlation coefficients at the best delay,
%	  with maxcorrthresh setting the cutoff [0.7 or 0.75 seem about right])
%	- If function "listen" is turned on, Figure 9 will be like Figure 7
%	  but just for a single analysis interval, and if 'showfull' = 1,
%	  Figure 10 will be like Figure 7 but will also show the location of
%	  the analysis window in the full signal.  (This can be slow to plot!)
%
%	HISTORY
%	- split from locatedriver.m into findfrogs.m, the main driver and
%	  initial correlation and delay processing, and extractfrogs.m,
%	  which locates recovers times and signals of individual frogs.
%
%	DEVELOPER NOTES
%	- correlating with smallest-power channel seems like a bad idea;
%	  maybe sectorize and use closest three mics for each sector?
%
% - 0.65 on maxcorrthresh seems to be getting close to the conservative limit;
%   0.7 seems very safe.  Down in the 0.63 range a few bad blocks are accepted.
%   Some down in the 0.4x range are OK, so if one is willing to accept some false
%   alarms, that could be used.
%
% - Should threshhold on the min delay error as well
%
% - regularization factor of 0.1 seemed to result in self- attenuation of less than 5%
%   Regularization of 0.05 gave self-attenuation of about 10%, but perhaps almost,
%   but not quite, as much attenuation of the interferers.  Using FFT length of 512
%   instead of 1024 (reg still 0.05) showed self-attenuation of about 15%, and overall
%   comparable but often worse at peaks attenuation of interferers; 1024, 0.1 probably
%   still preferred choice.  Haven't looked at the REALLY close ones, though, where
%   less regularization may help.
%
%  HYLA CINEREA NOTES
%
% - There seems to be some weird low-frequency baseline wobble in channels 7 and 8 after
%   filtering, which shouldn't be there ... check whether the transition bands are too small.
%
% - Filtering from 700 to 1500 Hz yields significantly higher xcorr values ...
%   maybe try limiting to this band for initial location estimate?
%
% - Channels 7 and 8 have significantly lower average xcorrs ... try leaving them out?
%
% - When applying the algorithm to loudspeaker data (max level about 0.7 (700 to 1500 filter))
%   maxcorrthresh = 0.65, all 8 channels, nails every pulse and only lets a few bad ones through
%   [700 1500] works better than [700 5500] here as well.
%
% - At least with the [700 1500] filter on the loudspeaker data, 0.1 sec corrblocksize
%   works better than 0.05 sec.  And 0.2 better still, in terms of reduced total delay
%   error
%
% - nitsgrad seems to need to be 100 or more occasionally
%
% - Is there an issue with including the reference channel (measured delay always zero)
%   in the total delay error metric? (I think not)
%
% - Looking at the results for the loudspeaker, which ought to be about as good as it gets,
%   the data just don't look very consistent.  Even the large-amplitude waveforms do not
%   look very much like delayed versions of each other, and the delays just don't line
%   up real well.  The maxcorrs are definitely higher for all channels when the signal
%   is present (all 8 above 0.5, often all above 0.6) as opposed to min of 0.2 or 0.3
%   when not.  I now suspect that multipath is killing us.  This will kill the location
%   estimate, but might not hurt extraction, as long as the apparent locations cluster.
%   Also looks somewhat more consistent for [700 5000] filter when signal is present,
%   but the threshholds drop a lot.  [2000 4000] is also more consistent than the full band,
%   although not as much so as the lower band.
%      Just looking at the spectrum, I see very different peaks and valleys for the
%   different channels above 1 kHz, which again suggests the multipath problem.
%      Ah, now [800 1050] gets interesting ... maxcorrthresh = 0.5, rejectlocthresh = 300,
%   we get 18 active frames, apparently in the right places.  Actually, looks like I could
%   go back to my old maxcorrthresh = 0.65 and get most of them.  The location seems off,
%   but the smear is only about a meter in each dimension.
%

%   SETUP

%   set parameters
close all; clc; clear;

showfull = 0;

csound = 343;       % m/sec
fs = 20000;         % Sampling rate

refchannel = 1;	% the channel to cross-correlate with for relative delay estimates
locchannelct = 5;
extractchannelct = 5;
rejectdelaysct = 2;
selectweight = 0.0;	% rates channels based on selectweight*powers + (1-selectweight)*maxcorrs

%startsamp = 0*fs+1;
startsamp = 40*fs+1;	% enter start_time_into_record_in_seconds*fs+1
endsamp = 55*fs;	% enter end_time_into_record_in_seconds*fs+1; if > 50, it crosses record boundary, and loadFileCrossFrames15Chan is used


%framesize = 0.02*fs;	% good numbers for cricket frogs
%corrblocksize = 0.03*fs;	% good numbers for cricket frogs
%framesize = 0.5*fs;	% good numbers for gulf coast toads
%corrblocksize = 0.5*fs;	% good numbers for toads
framesize = 0.04*fs;		% good numbers for Hyla cinerea?
corrblocksize = 0.12*fs;	% good numbers for Hyla cinerea?

%freqbands = [680 2300; 2700 4000];
%freqbands = [680 2300];	% Gulf Coast Toad band
%freqbands = [2700 4000];	% cricket frog band
freqbands = [700 5500];	% Hyla cinerea = Green Tree Frog
%freqbands = [2400 7000];	% Narrow Mouth Toad
%freqbands = [800 1000];
bandct = size(freqbands,1);

%miclocs = [0 0 0.11; 0 1 0.11; 1 1 0.11; 1 0 0.11];
%miclocs = [2 2 1.69; 2 5 1.74; 5 5 1.7625; 5 2 1.65];	% Cibolo omni locations
%miclocs = [0 0 1.7526; 0 0.1651 3.1242; -3.3020 0 1.7272; -3.3274 0.2540 3.0988; -3.3782 -1.9812 1.7399; -3.4155 -2.1255 3.1483; -0.0642 -1.9601 1.7960; -0.0822 -2.1152 3.1626];	% Day 2 BBSP locations
%miclocs = [ 3.0861  0  0.9398; 3.1115  0.1397  3.1064; 0  0  0.9906; 0  0.1334  3.1572; -3.0734  0  1.0033; -3.0680  0.1652  3.1248; -0.0313  -2.5252  0.9983; -0.0282  -2.6354  3.1754];	% Day 1 BBSP locations
%miclocs = [-0.5 0 3;.5 0 3;0 -.5 3;0 .5 3;0 0 2]; %simulation mic locations
%miclocs = [-1.5 -1.5 2.5; -1.5 1.5 2.5; 1.5 -1.5 2.5; 1.5 1.5 2.5] 3 by 3 square; centered at origin; 2.5 up
%miclocs = [0 .7071 2.63; .7071 0 2.63; 0 -.7071 2.63; -.7071 .0 2.63; 0 .7071 1.63] %1 by 1 square; centered at origin
%miclocs = [.79 3.94 2.80; 2.20 1.50 2.80; .79 2.53 2.80; 0.08 3.23 2.80; .79 3.94 1.80; 4.01 6.05 2.63; 4.72 5.35 2.63; 4.01 4.53 2.63; 3.30 5.35 2.63; 4.01 6.05 1.63; -2.51 7.41 2.91; -1.80 6.70 2.91; -2.51 5.99 2.91; -3.22 6.70 2.91; -2.51 7.41 1.91] %2009 Day 2 from Russell; erroneous!
miclocs = [    0.7900    4.0371    2.8034
    1.4971    3.3300    2.8034
    0.7900    2.6229    2.8034
    0.0829    3.3300    2.8034
    0.7900    4.0371    1.8034
    4.0100    7.0371    2.6256
    4.7171    6.3300    2.6256
    4.0100    5.6229    2.6256
    3.3029    6.3300    2.6256
    4.0100    7.0371    1.6256
   -2.5000    9.0471    2.9050
   -1.7929    8.3400    2.9050
   -2.5000    7.6329    2.9050
   -3.2071    8.3400    2.9050
   -2.5000    9.0471    1.9050]; % Night 2, BBSP 2009, calculated by arraylocations.m
%  On 6/18-19/09, Pole of Array 1 at (0.79,3.23,1.80) (z = vertical height of lower element)
% Pole of Array 2 at (4.01,5.35,1.63)
% Pole of Array 3 at (-2.51,6.70,1.91)
micct = size(miclocs,1);
dimspace = size(miclocs,2);
%miclocs = [0 0 0.11; 0 1 0.11; 1 1 0.11];

fftlen = 1024; % length of FFT to use in FMV adaptive beamformer for recovering individual source; should be several times longer than maxdelay.  2048 does marginally better than 1024 for frogs 4+5 under the others, maybe because corrlen is longer.  Seems to be less self-cancellation w 2048?
fftblockstep = fftlen/8; % stepsize for steering vector estimates
fmindex = ceil(fftlen*freqbands(:,1)/fs)
fmaxindex = floor(fftlen*freqbands(:,2)/fs)

%  Set grid to search for sources

xmin = -4;
xmax = 12;
xinc = 0.1;
xct = floor((xmax-xmin)/xinc)+1;
ymin = -0;
ymax = 14;
yinc = 0.1;
yct = floor((ymax-ymin)/yinc)+1;
zmin = -0.5;
zmax = 1;

maxcorrthresh = 0.3;	% 0.65 is max aggressive; start to get a lot of total outliers below that.  0.75 leaves very few outliers, but keeps most good ones
powersratiothresh = 0.05;
	% 0.02 is max aggressive; below that, toads start to bleed through into cricket frog band.  0.08 makes for a very nice picture of cricket frogs
rejectlocthresh = 10;	% if mean errors in delay exceeds this threshhold, reject location estimate.  For toads with 4 mics, 10-20 seemed to be good; things about 50 bad

coherence_test_flag = 0;	% Check whether sum of largest eigenvalues over band exceeds threshold
coherencethresh = 0.5;

fftlencoherence = fftlen;
coherence_step = fftlencoherence/4;
fmindexcoh = ceil(fftlencoherence*freqbands(:,1)/fs);
fmaxindexcoh = floor(fftlencoherence*freqbands(:,2)/fs);

%w = ones(micct,fftlencoherence);  % make a rectangular window
%wb = hann(fftlencoherence);%hamming(winlen);
%w = zeros(fftlencoherence,micct);	%  Make the STFT window
%for ii=1:channelct,
%  w(:,ii) = [zeros(1,(binct-winlen)/2) wb' zeros(1,(binct-winlen)/2)]'/max(wb);
%  w(:,ii) = hann(binct);
%end

locbygradflag = 1;
locbypowerflag = 0;
gradcompstep = 0.0001;
mustart = 1;
mumin = 0.000001;
mumax = 10;
maxstepsq = 2.0^2;
muup = 2;
mudown = 1/3;
nitsgrad = 125;

%
%  Create gain and delay templates for each location in search grid to compare with actual data
%
if (locbygradflag < 1),
    locgrid = zeros(xct*yct,dimspace);
    locct = xct*yct;
    
    iiii = 0;
    for xloc=xmin:xinc:xmax,
        for yloc=ymin:yinc:ymax,
            iiii = iiii + 1;
            locgrid(iiii,:) = [xloc yloc 0];
        end
    end
    [gains,delays] = makedd(miclocs,locgrid,fs);
end

%
%  LOAD THE DATA
%
%load 070321e_002.mat	% Cibolo load

%%% Load 2008 BBSP data
%cd BBSP  % Load for 2008 BBSP
%loadFileBBSP
%cd ..
%x = micarray(startsamp:endsamp,1:micct)';
%%%
%
%Load 12 Channel Data
%loadFile12Chan
%x = micarray(startsamp:endsamp,9:12)';
%
%%% Load 15 Channel Data
cd BBSP090619
if (endsamp <= 50*fs),
   loadFile15Chan
else
   loadFileCrossFrames15Chan
   clear micarray1
end
cd ..
%%% End of 15-channel load code

%whichmics = [5 9 10 11 15]  % Good set for frog E? Seems high in the air
%whichmics = [1 5 10 13 15]  % Good for frogs A-D, 0190/0191
%whichmics = [5 6 7 8 15]  %  Best for frog F, 0195/0196
%whichmics = [1 2 4 6 8]  % high mics + 1 for BBSP 2008 data
%whichmics = [1 3 5 7 6] % low mics + 1 for BBSP 2008 data
%x = micarray(startsamp:endsamp,whichmics)';
%miclocs = miclocs(whichmics,:);

%%% Code for 15-channel, three-array data from BBSP 2009
whicharray = input('Recover with 1 = Box 1 array; 2 = Box 2; 3 = Box 3?  ')
x = micarray(startsamp:endsamp,5*(whicharray-1)+1:5*whicharray)';
miclocs = miclocs(5*(whicharray-1)+1:5*whicharray,:)
clear micarray;
pack

figure(1);
scatter3(miclocs(:,1), miclocs(:,2), miclocs(:,3));
%%%

micct = size(miclocs,1);

%
%  Compute distances from estimated locations to set maximum delay
%
distfromlocs = zeros(micct,micct);
for ii=1:micct,
   for iii=1:micct,
      distfromlocs(ii,iii) = sqrt(sum((miclocs(ii,:)-miclocs(iii,:)).^2));
   end
end
maxdelay = ceil(max(max(distfromlocs))*fs/csound + 10)%+ 30)
mindelay = - maxdelay;

%
%Load simfrog data
%x = signal; 
xlen = max(size(x));

tottimes = [1:xlen]/fs;

plotccout(x,zeros(size(tottimes)),tottimes,10)	% Plot the raw data in figure 10

%
%  COMPUTE THE DELAYS AND POWERS
%
framect = floor((xlen-corrblocksize-maxdelay+mindelay)/framesize) + 1;
measdelays = zeros(bandct,micct,micct,framect);
maxcorrs = zeros(bandct,micct,micct,framect);
powers = zeros(bandct,micct,framect);
coherence = ones(bandct,framect);
mindelaydiffs = zeros(bandct,framect);
mindelaylocs = zeros(bandct,framect);
mindelays = zeros(bandct,micct,framect);
minlocs = zeros(bandct,framect,dimspace);
refchannels = zeros(bandct,framect);
channelselectmetrics = zeros(micct,framect);

w = ones(micct,fftlencoherence);  % make a rectangular window for coherence test
%wb = hann(fftlencoherence);%hamming(winlen);
%w = zeros(fftlencoherence,micct);	%  Make the STFT window
%for ii=1:channelct,
%  w(:,ii) = [zeros(1,(binct-winlen)/2) wb' zeros(1,(binct-winlen)/2)]'/max(wb);
%  w(:,ii) = hann(binct);
%end

%%
for ii=1:bandct, %ii is the band, always 1 for now
    xfilt = filterit(x,freqbands(ii,1),freqbands(ii,2),fs);	% filter to extract band of interest
    plotccout(xfilt,zeros(size(tottimes)),tottimes,10+ii)	% Plot the filtered data in figure 11
    %a = xfilt(3,:);
    %save "call.mat" a %save data for simfrog
    %pause
    iii = 0;
    % loop through frames
    for iiii=1-mindelay:framesize:xlen-corrblocksize-maxdelay,
        iii = iii + 1 % iii: iframe
        % compute energy of frame in each channel
        for jjj=1:micct,	
            powers(ii,jjj,iii) = xfilt(jjj,iiii:iiii+corrblocksize-1)* xfilt(jjj,iiii:iiii+corrblocksize-1)'/corrblocksize;
            maxcorrs(ii,jjj,jjj,iii) = 1.0;	% correlation with self always 1
        end
        [maxframepower,refchannel] = max(squeeze(powers(ii,:,iii)))
        refchannels(iii) = refchannel;% get the channel with highest power
        % compute normalize cross correlation between pair of mic from
        % mindelay to maxdely
        for jj=1:micct-1,
            for jjj=jj+1:micct,
                % ****tNOTE: powers value are overwritely updated******
                [measdelays(ii,jj,jjj,iii),maxcorrs(ii,jj,jjj,iii),powers(ii,jj,iii),powers(ii,jjj,iii)] = xcorrdat(xfilt(jj,:),xfilt(jjj,:),iiii,corrblocksize,mindelay,maxdelay);
                measdelays(ii,jjj,jj,iii) = - measdelays(ii,jj,jjj,iii);	% same result so symmetrize
                maxcorrs(ii,jjj,jj,iii) = maxcorrs(ii,jj,jjj,iii);	% symmetric
            end
        end
        %  Compute max-eigenvalue-based coherence test
        if (coherence_test_flag == 1),  
            % Compute short-time-frequency correlation matrices
            Rstcoh = zeros(fftlencoherence,micct,micct);
            for iiiii=iiii:coherence_step:iiii+corrblocksize-fftlencoherence+1,
                fftbuf = fft(w.*x(:,iiiii:iiiii+fftlencoherence-1),[],1);
                for iiiiii=fmindexcoh(ii):fmaxindexcoh(ii),
                    Rstcoh(iiiiii,:,:) = squeeze(Rstcoh(iiiiii,:,:)) + fftbuf(:,iiiiii)*fftbuf(:,iiiiii)';	% recursively update running correlations in each bin
                end
            end
            
            maxeigsum = 0;
            totaleigsum = 0;
            for iiiiii=fmindexcoh(ii):fmaxindexcoh(ii),	% Compute correlation metric over band
                Rxcorruse = squeeze(Rstcoh(iiiiii,:,:));
                D = real(eig(Rxcorruse));
                maxeigsum = maxeigsum + max(D);
                totaleigsum = totaleigsum + sum(D);
            end
            coherence(ii,iii) = maxeigsum/totaleigsum;
        end
        channelselectmetrics(:,iii) = (selectweight*squeeze(powers(ii,:,iii))/maxframepower + (1-selectweight)*squeeze(maxcorrs(ii,refchannels(iii),refchannels(iii),iii)))';
    end % end loop through frames
end
%%
%
%  GENERATE THE LOCATION MAPS
%

%realdelays = [204960 204880 204960 205010]';	% from 070321e_002  2.5 4.5 toad
%realgains = [0.62 0.70 0.52 0.40]';	% from 070321e_002

%realdelays = [77586 77472 77620 77663]';	% from 070321e_002   loser toad
%realgains = [0.36 0.52 0.35 0.31]';	% from 070321e_002

%
%  FIND AND DISPLAY THE LOCATIONS
%
showframesflag = 0;
if (locbygradflag == 0), showframesflag = input('Show location map for each frame?  1=yes, 0=no  '); end

%%
for ii=1:bandct, %ii is the band, always 1 for now
    %%
    ii = 1;
    maxpowers = max(max(squeeze(powers(ii,:,:))));
    for iiii=1:framect, %iiii is the frame
        
        iiii
        realdelays = squeeze(measdelays(ii,:,:,iiii));%';
        %    realdelays = realdelays-ones(micct,1)*mean(mean(realdelays)); don't remove mean for purely relative measurements!
        realgains = sqrt(squeeze(powers(ii,:,iiii)))';
        
        
        %    gaincorr = (gains'*(realgains/sqrt(realgains'*realgains)))';
        %    [maxgaincorr,maxgainindex] = max(gaincorr)
        if ( (showframesflag == 0) & (locbygradflag == 1) & (locbypowerflag == 0) ),
            %
            %  LOCATE VIA GRADIENT DESCENT
            %
            [junk,micsort] = sort(channelselectmetrics(:,iiii));
            miclocsbest = miclocs(micsort(micct-locchannelct+1:micct),:);
            startloc = mean(miclocsbest);	% start in center of working microphone array
            %     startloc(3) = 0;	% start at water surface
            realdelays = realdelays(micsort(micct-locchannelct+1:micct),micsort(micct-locchannelct+1:micct));	% just include the "best" set of microphones
            
            maskmat = ones(size(realdelays));	% matrix for masking out rejected delays
            for i7=0:rejectdelaysct,	% loop for leaving out bad relative delays
                locnow = startloc;
                mu = mustart;
                for iiiii=1:nitsgrad,
                    [junk,delaysnow] = makedd(miclocsbest,locnow,fs);	% find delays for current best location estimate
                    delaysrelativemat = ones(locchannelct,1)*delaysnow' - delaysnow*ones(1,locchannelct);
                    %              realdelays
                    %              pause
                    meas0 = mean(mean((maskmat.*(delaysrelativemat - realdelays)).^2));	% total squared delay error is metric to minimize
                    for iiiiii=1:dimspace,	% compute numerical gradient
                        locnowtemp = locnow;
                        locnowtemp(iiiiii) = locnowtemp(iiiiii) + gradcompstep;
                        [junk,delaysnowtemp] = makedd(miclocsbest,locnowtemp,fs);
                        delaysrelativemat = ones(locchannelct,1)*delaysnowtemp' - delaysnowtemp*ones(1,locchannelct);
                        meas = mean(mean((maskmat.*(delaysrelativemat - realdelays)).^2));	% total squared delay error is metric to minimize
                        gradxyz(iiiiii) = (meas - meas0)/gradcompstep;
                    end
                    if (mu*gradxyz*gradxyz' > maxstepsq),	% limit stepsize
                        mu = min([maxstepsq/sqrt(gradxyz*gradxyz') mu mumax]);
                        %         info = 'exceeded max stepsize'
                    end
                    % gradient descend update
                    locnowtemp = locnow - mu*gradxyz;
                    if (locnowtemp(1) < xmin), locnowtemp(1) = xmin; end	% enforce outer region boundaries
                    if (locnowtemp(1) > xmax), locnowtemp(1) = xmax; end
                    if (locnowtemp(2) < ymin), locnowtemp(2) = ymin; end
                    if (locnowtemp(2) > ymax), locnowtemp(2) = ymax; end
                    if (locnowtemp(3) < zmin), locnowtemp(3) = zmin; end
                    if (locnowtemp(3) > zmax), locnowtemp(3) = zmax; end
                    % locnowtemp(3) = 0;	% projection for frogs only on flat surface
                    [junk,delaysnowtemp] = makedd(miclocsbest,locnowtemp,fs);
                    delaysrelativemat = ones(locchannelct,1)*delaysnowtemp' - delaysnowtemp*ones(1,locchannelct);
                    meas = mean(mean((maskmat.*(delaysrelativemat - realdelays)).^2));	% total squared delay error is metric to minimize
                    if (meas < meas0),
                        locnow = locnowtemp;
                        meas0 = meas;
                        mu = min(mu*muup,mumax);
                        %info = 'accept new location; increase mu'
                    else
                        mu = mu*mudown/muup;
                        %info = 'reject new location; lower mu'
                    end
                    %if (mu < mumin), break; end
                    %iiiii
                    %locnow
                    %gradxyz
                    %muF
                    %meas
                    %meas0
                    %pause
                end % end loop for gradient descend
                if (i7 < rejectdelaysct),
                    [temp,index1] = max(abs(maskmat.*(delaysrelativemat - realdelays)));	% find element of largest deviation
                    [junk,index2] = max(abs(temp));
                    maskmat(index1(index2),index2) = 0;	% mask out element of largest deviation
                    maskmat(index2,index1(index2)) = 0;
                    %              maskmat
                end
            end
            
            minlocs(ii,iiii,:) = locnow;		% save minimizing location for this frame
            [junk,delaysnow] = makedd(miclocsbest,locnow,fs);	% find delays for current best location estimate
            delaysrelativemat = ones(locchannelct,1)*delaysnow' - delaysnow*ones(1,locchannelct); %BUG: delaysnow instead of delaysnowtemp
            mindelaydiffs(ii,iiii) = mean(mean((abs(maskmat.*(delaysrelativemat - realdelays)))));
            mindelays(ii,:,iiii) = delaysnow;
            %         realdelays
            %         maskmat
            %         pause
        elseif ( (showframesflag == 0) & (locbygradflag == 1) & (locbypowerflag == 1) ),
            [junk,micsort] = sort(channelselectmetrics(:,iiii));
            miclocsbest = miclocs(micsort(micct-locchannelct+1:micct),:);
            startloc = mean(miclocsbest);	% start in center of working microphone array
            realgains = realgains(micsort(micct-locchannelct+1:micct,:));
            realgainsinv = ones(size(realgains))./realgains;
            startloc(3) = 0.0;	%  all frogs are on the ground at z=0 elevation at Cibolo
            locnow = startloc;%[ 1 1 0 ]
            mu = mustart;
            for iiiii=1:nitsgrad,
                [gainsnow,delaysnow] = makedd(miclocsbest,locnow,fs);	% find delays and gains for current best location estimate
                gainsnowinv = ones(size(gainsnow))./gainsnow;
                alpha = (realgainsinv'*gainsnowinv)/(gainsnowinv'*gainsnowinv);
                meas0 = (realgainsinv - alpha*gainsnowinv)'*(realgainsinv - alpha*gainsnowinv);	% squared inverted gains error is metric to minimize
                for iiiiii=1:dimspace,	% compute numerical gradient
                    locnowtemp = locnow;
                    locnowtemp(iiiiii) = locnowtemp(iiiiii) + gradcompstep;
                    [gainsnowtemp,delaysnowtemp] = makedd(miclocsbest,locnowtemp,fs);
                    gainsnowtempinv = ones(size(gainsnowtemp))./gainsnowtemp;
                    alpha = (realgainsinv'*gainsnowtempinv)/(gainsnowtempinv'*gainsnowtempinv);
                    meas = (alpha*gainsnowtempinv - realgainsinv)'*(alpha*gainsnowtempinv - realgainsinv);
                    gradxyz(iiiiii) = (meas - meas0)/gradcompstep;
                end
                if (mu*gradxyz*gradxyz' > maxstepsq),	% limit stepsize
                    mu = min([maxstepsq/sqrt(gradxyz*gradxyz') mu mumax]);
                    %info = 'exceeded max stepsize'
                end
                locnowtemp = locnow - mu*gradxyz;
                if (locnowtemp(1) < xmin), locnowtemp(1) = xmin; end	% enforce outer region boundaries
                if (locnowtemp(1) > xmax), locnowtemp(1) = xmax; end
                if (locnowtemp(2) < ymin), locnowtemp(2) = ymin; end
                if (locnowtemp(2) > ymax), locnowtemp(2) = ymax; end
                if (locnowtemp(3) < zmin), locnowtemp(3) = zmin; end
                if (locnowtemp(3) > zmax), locnowtemp(3) = zmax; end
                [gainsnowtemp,delaysnowtemp] = makedd(miclocsbest,locnowtemp,fs);
                gainsnowtempinv = ones(size(gainsnowtemp))./gainsnowtemp;
                alpha = (realgainsinv'*gainsnowtempinv)/(gainsnowtempinv'*gainsnowtempinv);
                meas = (alpha*gainsnowtempinv - realgainsinv)'*(alpha*gainsnowtempinv - realgainsinv);
                if (meas < meas0),
                    locnow = locnowtemp;
                    meas0 = meas;
                    mu = min(mu*muup,mumax);
                    % info = 'accept new location; increase mu'
                else
                    mu = mu*mudown/muup;
                    % info = 'reject new location; lower mu'
                end
            end
            minlocs(ii,iiii,:) = locnow;		% save minimizing location for this frame
            [gainsnow,delaysnow] = makedd(miclocsbest,locnow,fs);	% find delays for current best location estimates
            alpha = (gainsnow'*delaysnow)/(delaysnow'*delaysnow);
            mingaindiffs(ii,iiii) = sum(abs(gainsnow*alpha - realgains));
            [gainsnow,delaysnow] = makedd(miclocs,locnow,fs);	% find delays for current best location estimate
            mindelays(ii,:,iiii) = delaysnow;
        else
            %
            %  LOCATE VIA GRID SEARCH
            %
            delaydiffs = zeros(1,locct);
            for iii=1:locct,
                delaydiffs(iii) = sum(abs(delays(:,iii) - realdelays));
            end
            
            time = (iiii-1)*framesize/fs
            if (showframesflag == 1), figure(1); end
            %    imagesc(flipud(reshape(gaincorr,xct,yct)))
            max(squeeze(powers(ii,:,iiii)));
            if (max(squeeze(powers(ii,:,iiii))) < powersratiothresh*maxpowers),
                temp = zeros(size(delaydiffs));
                temp(1) = 1;
                if (showframesflag == 1), imagesc(flipud(reshape(temp,yct,xct))); end
            elseif (min(min(squeeze(maxcorrs(ii,:,:,iiii)))) < maxcorrthresh),
                if (showframesflag ==1), imagesc(flipud(reshape(ones(size(delaydiffs)),yct,xct))); end
            else
                if (showframegainsnowsflag == 1), imagesc(flipud(reshape(delaydiffs,yct,xct))); end
                [mindelaydiffs(ii,iiii),mindelaylocs(ii,iiii)] = min(delaydiffs);
                minlocs(ii,iiii,:) = locgrid(mindelaylocs(ii,iiii),:);
            end
            
            if (showframesflag == 1), figure(2); end
            %    imagesc(flipud(reshape(delaydiffs,xct,yct)))
            if (showframesflag == 1),
                imagesc(flipud(reshape(ones(size(delaydiffs))./delaydiffs,yct,xct)))
                squeeze(maxcorrs(ii,:,:,iiii))
                squeeze(powers(ii,:,iiii))
                %    figure(3)
                %     listen(xfilt',(1+(iiii-1)*framesize-mindelay)/fs,((iiii-1)*framesize+corrblocksize-mindelay)/fs,1,4,fs,showfull)
                pause
            end
        end
    end % end frames loops framect
    
    %% Make cumulative plot of all locations
    command = 1;
    while (command ~= 10)
        figure(20)
        clf
        hold on
        plot(miclocs(:,1),miclocs(:,2),'ro')
        xdithermag = 0.1;
        ydithermag = 0.1;
        if (locbygradflag == 0),
            xdithermag = xinc;
            ydithermag = yinc;
        end
        activeframect = 0;
        for iii=1:framect,
            xdither = xdithermag*(rand(1)-0.5);
            ydither = ydithermag*(rand(1)-0.5);
            [junk,micsort] = sort(channelselectmetrics(:,iii));
            micsortbest = micsort(micct-locchannelct+1:micct);
            if ((max(squeeze(powers(ii,micsortbest,iii))) > powersratiothresh*maxpowers) &...
                    (min(squeeze(maxcorrs(ii,micsortbest,:,iii))) > maxcorrthresh) & ...
                    (mindelaydiffs(ii,iii) < rejectlocthresh) & (coherence(ii,iii) > coherencethresh) ), %(mindelaylocs(ii,iii) > 0),
                plot(minlocs(ii,iii,1)+xdither,minlocs(ii,iii,2)+ydither,'+')
                activeframect = activeframect + 1;
            end
        end
        hold off
        grid on
        
        activeframect
        maxcorrthresh
        powersratiothresh
        rejectlocthresh
        coherencethresh
        command = input('COMMANDS: \n 1 Decrease Maximum Correlation Threshold \n 2 Increase Maximum Correlation Threshold \n 3 Decrease Power Ratio Threshold \n 4 Increase Power Ratio Threshold \n 5 Decrease Location Threshold \n 6 Increase Location Threshold \n 7 Decrease Coherence Threshold \n 8 Increase Coherence Threshold \n 9 Enter Threshold Values \n 10 Done \n>>  ');
        if command == 1
            maxcorrthresh = maxcorrthresh^(10/9)
        elseif command == 2
            maxcorrthresh = maxcorrthresh^.9
        elseif command == 3
            powersratiothresh = powersratiothresh^1.25
        elseif command == 4
            powersratiothresh = powersratiothresh^0.8
        elseif command == 5
            rejectlocthresh = rejectlocthresh*2/3
        elseif command == 6
            rejectlocthresh = rejectlocthresh*3/2
        elseif command == 7
            coherencethresh = coherencethresh^(10/9)
        elseif command == 8
            coherencethresh = coherencethresh^0.9
        elseif command == 9
            a = maxcorrthresh
            maxcorrthresh = input('Maximum Correlation Threshold: ');
            if maxcorrthresh == -1, maxcorrthresh = a; end
            a = powersratiothresh
            powersratio = input('Power Ratio Threshold: ');
            if powersratiothresh == -1, powersratiothresh = a; end
            a = rejectlocthresh
            rejectlocthresh = input('Reject Location Threshold: ');
            if rejectlocthresh == -1, rejectlocthresh = a; end
            a = coherencethresh
            coherencethresh = input('Coherence Threshold: ');
            if coherehcethresh == -1, coherencethresh = a; end
        end
    end
    %framect
    %activeframect
    %rejectlocthresh
    %squeeze(maxcorrs)
    %squeeze(mindelaydiffs)
    %delaydiffsrel = squeeze(mindelays) - squeeze(measdelays);
    %delaydiffsrel = delaydiffsrel - ones(micct,1)*delaydiffsrel(refchannel,:)
    
    activeframect
    %
    % EXTRACT INDIVIDUAL SOURCES FROM SIGNAL
    %
%     extractfrogs3D
    
end


%   DONE!
