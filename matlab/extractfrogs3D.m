%
%		extractfrogs3D.m
%		Douglas L. Jones
%		University of Illinois
%		September 10, 2008
%
%	extractfrogs3D.m: script that takes frame-by-frame cross-correlation
%		and location data computed in findfrogs.m,
%		shows times and locations of solo sources, estimates steering vectors
%		and recovers individual sources using FMV
%
%	USER NOTES
%	- Figure 20 plots the locations of the source from each reliable frame
%	- Figures 21 - 20+whichsource plot the locations and the bounding box
%	  the user selects for a specific source
%	- Figures 31-30+whichsource plots the four omni signal channels and
%	  the time-locations at which the selected source occurs
%	- Figure 41 - 40+whichsource displays that source as recovered by FMV,
%	  either from the omni array or from the XYZ0 colocated array
%
%	HISTORY
%	- extractfrogs.m was broken out from locatedriver, so that it can
%	  be run multiple times without reprocessing from the start
%  - extractfrogs3D.m extends extractfrogs to be compatible with findfrogs3D.m
%    and to add some new features
%
%	DEVELOPER NOTES
%

%
% EXTRACT INDIVIDUAL SOURCES FROM SIGNAL
%
sourceextractct = input('How many sources would you like to recover?  ')

thissource = zeros(sourceextractct,xlen);
whentimes = zeros(sourceextractct,xlen);
lowcorners = zeros(sourceextractct,2);
highcorners = zeros(sourceextractct,2);
boundingbox = zeros(5,2);
extractxyzflag = zeros(sourceextractct,1);

meanlocs = zeros(sourceextractct,3);
stdlocs = zeros(sourceextractct,1);

for whichsource=1:sourceextractct,
    figure(20+whichsource)
    if (locbygradflag == 0),
      xdithermag = xinc;
      ydithermag = yinc;
    end
    clf
    hold on
    plot(miclocs(:,1),miclocs(:,2),'ro')
    activeframect = 0;
    for iii=1:framect,
      xdither = xdithermag*(rand(1)-0.5);
      ydither = ydithermag*(rand(1)-0.5);
    
      [ans,micsort] = sort(channelselectmetrics(:,iii));
      micsortbest = micsort(micct-locchannelct+1:micct);
      if ((max(squeeze(powers(ii,micsortbest,iii))) > powersratiothresh*maxpowers) & (min(squeeze(maxcorrs(ii,micsortbest,:,iii))) > maxcorrthresh) & (mindelaydiffs(ii,iii) < rejectlocthresh) & (coherence(ii,iii) > coherencethresh) ), %(mindelaylocs(ii,iii) > 0)
        plot(minlocs(ii,iii,1)+xdither,minlocs(ii,iii,2)+ydither,'+')
        activeframect = activeframect + 1;
      end
    end
    hold off
    grid on
    
    parameter_change_flag = 0;
    %
    %  CAN'T USE THIS FEATURE UNTIL YOU SAVE INDIVIDUAL PARAMETERS FOR EACH FROG AND USE THEM WHEN FORMING INDICATOR FUNCTIONS FOR STEERING VECTOR ESTIMATION!!
    %
    %parameter_change_flag = input('Would you like to change thresholds for source location clustering?  0=no, 1=yes  ')
    if ( parameter_change_flag >= 1 ),
       % Make cumulative plot of all locations
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
			if ((max(squeeze(powers(ii,micsortbest,iii))) > powersratiothresh*maxpowers) & (min(squeeze(maxcorrs(ii,micsortbest,:,iii))) > maxcorrthresh) & (mindelaydiffs(ii,iii) < rejectlocthresh) & (coherence(ii,iii) > coherencethresh) ), %(mindelaylocs(ii,iii) > 0),
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
   end
    
 
    %plotwhere(minlocs,powers,maxcorrs,framect,ii,maxpowers,powersratiothresh,maxcorrthresh,xdithermag,ydithermag,20+whichsource)

    boundflag = 0;
    while (boundflag ~= 1)
      whichsource
      lowcorners(whichsource,1) = input('Enter lower X corner coordinate of rectangle bounding this source  ')
      lowcorners(whichsource,2) = input('Enter lower Y corner coordinate of rectangle bounding this source  ')
      highcorners(whichsource,1) = input('Enter upper X corner coordinate of rectangle bounding this source  ')
      highcorners(whichsource,2) = input('Enter upper Y corner coordinate of rectangle bounding this source  ')
      %
      %  Plot the bounding box for confirmation
      %
      boundingbox(1,:) = lowcorners(whichsource,:);
      boundingbox(2,:) = [lowcorners(whichsource,1) highcorners(whichsource,2)];
      boundingbox(3,:) = highcorners(whichsource,:);
      boundingbox(4,:) = [highcorners(whichsource,1) lowcorners(whichsource,2)];
      boundingbox(5,:) = lowcorners(whichsource,:);
      clf
      hold on
      plot(miclocs(:,1),miclocs(:,2),'ro')
      
      activeframect = 0;
      for iii=1:framect,
        xdither = xdithermag*(rand(1)-0.5);
        ydither = ydithermag*(rand(1)-0.5);
        
        [ans,micsort] = sort(channelselectmetrics(:,iii));
        micsortbest = micsort(micct-locchannelct+1:micct);
	     if ((max(squeeze(powers(ii,micsortbest,iii))) > powersratiothresh*maxpowers) & (min(squeeze(maxcorrs(ii,micsortbest,:,iii))) > maxcorrthresh) & (mindelaydiffs(ii,iii) < rejectlocthresh) & (coherence(ii,iii) > coherencethresh) ), %(mindelaylocs(ii,iii) > 0),
          plot(minlocs(ii,iii,1)+xdither,minlocs(ii,iii,2)+ydither,'+')
          if ((minlocs(ii,iii,1) >= lowcorners(whichsource,1)) & (minlocs(ii,iii,2) >= lowcorners(whichsource,2)) & (minlocs(ii,iii,1) <= highcorners(whichsource,1)) & (minlocs(ii,iii,2) <= highcorners(whichsource,2))), activeframect = activeframect + 1; end
        end
      end
      plot(boundingbox(:,1)',boundingbox(:,2)','k')
      hold off
      grid on
      activeframect
      boundflag = input('Is this OK?  1 = yes, 0 = no  ');
    end
    %
    %  Show locations of these frames in signal
    %
    activeframect = 0;
    meanloc = zeros(1,dimspace);		% mean location of the selected source
    secondmomentloc = zeros(1,dimspace);
    for iii=1:framect,
      [ans,micsort] = sort(channelselectmetrics(:,iii));
      micsortbest = micsort(micct-locchannelct+1:micct);
      if ((max(squeeze(powers(ii,micsortbest,iii))) > powersratiothresh*maxpowers) & (min(squeeze(maxcorrs(ii,micsortbest,:,iii))) > maxcorrthresh) & (mindelaydiffs(ii,iii) < rejectlocthresh) & (coherence(ii,iii) > coherencethresh) ), %(mindelaylocs(ii,iii) > 0),
        if ((minlocs(ii,iii,1) >= lowcorners(whichsource,1)) & (minlocs(ii,iii,2) >= lowcorners(whichsource,2)) & (minlocs(ii,iii,1) <= highcorners(whichsource,1)) & (minlocs(ii,iii,2) <= highcorners(whichsource,2))),
          activeframect = activeframect + 1;
          meanloc = meanloc + squeeze(minlocs(ii,iii,:))';
          secondmomentloc = secondmomentloc + (squeeze(minlocs(ii,iii,:))').^2;
          whentimes(whichsource,(iii-1)*framesize-mindelay+1:(iii-1)*framesize-mindelay+corrblocksize) = ones(1,corrblocksize);
        end
      end
    end
    activeframect
    meanloc = meanloc/activeframect
    varloc = secondmomentloc/(activeframect) -meanloc.*meanloc;
    stdsxyz = sqrt(varloc)
    stdloc = sqrt(sum(varloc))
    meanlocs(whichsource,:) = meanloc;
    stdlocs(whichsource) = sqrt(sum(varloc));
%%%
%%% Code for saving cluster point locations
%    clusterlocs = zeros(activeframect,3);
%    activeframect = 0;
%    for iii=1:framect,
%      [ans,micsort] = sort(channelselectmetrics(:,iii));
%      micsortbest = micsort(micct-locchannelct+1:micct);
%      if ((max(squeeze(powers(ii,micsortbest,iii))) > powersratiothresh*maxpowers) & (min(squeeze(maxcorrs(ii,micsortbest,:,iii))) > maxcorrthresh) & (mindelaydiffs(ii,iii) < rejectlocthresh) & (coherence(ii,iii) > coherencethresh) ), %(mindelaylocs(ii,iii) > 0),
%         if ((minlocs(ii,iii,1) >= lowcorners(whichsource,1)) & (minlocs(ii,iii,2) >= lowcorners(whichsource,2)) & (minlocs(ii,iii,1) <= highcorners(whichsource,1)) & (minlocs(ii,iii,2) <= highcorners(whichsource,2))),
%           activeframect = activeframect + 1;
%           clusterlocs(activeframect,:) = minlocs(ii,iii,:);
%         end
%      end
%    end
%    if whichsource == 1
%       save sourcepoints1.mat clusterlocs whichmics meanloc stdsxyz stdloc
%    elseif whichsource == 2
%       save sourcepoints2.mat clusterlocs whichmics meanloc stdsxyz stdloc
%    elseif whichsource == 3
%       save sourcepoints3.mat clusterlocs whichmics meanloc stdsxyz stdloc
%    elseif whichsource == 4
%       save sourcepoints4.mat clusterlocs whichmics meanloc stdsxyz stdloc
%    elseif whichsource == 5
%       save sourcepoints5.mat clusterlocs whichmics meanloc stdsxyz stdloc
%    elseif whichsource == 6
%       save sourcepoints6.mat clusterlocs whichmics meanloc stdsxyz stdloc
%    end
%%%
%%%

      
    figure(30+whichsource)
    plotccout(xfilt,zeros(size(tottimes)),tottimes,30+whichsource)	% Plot the filtered data in figure 11
%    hold off
%    plot(tottimes,xfilt(1,:),'c')
    hold on
%    plot(tottimes,xfilt(2,:),'b')
%    plot(tottimes,xfilt(3,:),'g')
%    plot(tottimes,xfilt(4,:),'r')
    plot(tottimes,whentimes(whichsource,:)*1.1*max(max(abs(xfilt(:,:)))),'k')
    hold off
%    grid on
    %
    %  Find out which sources to extract
    %
    %extractxyzflag(whichsource) = input('Recover this source with  1=XYZO array, or 2=omni array, or 0=not at all?  ')
    extractxyzflag(whichsource) = input('Recover this source? 2=yes, 0=no  ')
    if (extractxyzflag(whichsource) > 0), extractxyzflag(whichsource) = 2; end
  end
  
%
%  EXTRACT SOURCES OF INTEREST
%
for whichsource=1:sourceextractct,
  whichsource
  if (extractxyzflag(whichsource) > 0),

    %
    %  Estimate steering vectors for this source
    %
    
    %  Find relative delays to microphones for this source to nearest sample
    [ans,meandelay] = makedd(miclocs,meanlocs(whichsource,:),fs)
    [ans,meandelayindex] = min(meandelay)
    refchannelindex = meandelayindex(1)
    meandelay = meandelay - meandelay(refchannelindex)*ones(size(meandelay))
    meandelay = round(meandelay)
    
    %  Apply temporal window (whentimes) to produce signals when that source only is active
    %  and delay all channels to line up to nearest sample
    if (extractxyzflag(whichsource) == 1)	% Recover using the colocated directional array
      refchannelindex = 4;
      xtemp = filterit(micarray(startsamp:endsamp,5:8)',freqbands(ii,1),freqbands(ii,2),fs);
      for iii=1:4,
	xtemp(iii,:) = xtemp(iii,:).*whentimes(whichsource,:);
      end
      if (whichsource == 1)
        figure(8)
        plot(xtemp')
        title('XYZO signals')
        grid on
      end
    else		% Use the spatially separated array to recover this frog
      xtemp = zeros(size(xfilt));
      for iii=1:micct,
        xtemp(iii,1:xlen-meandelay(iii)) = xfilt(iii,meandelay(iii)+1:xlen).*whentimes(whichsource,1:xlen-meandelay(iii));
      end
    end
    
    %  In each good frame, FFT data, cross-correlate in each frequency bin, and sum into correlation estimate
    Rxcorr = zeros(fftlen,micct,micct);
    channelselectmean = zeros(size(channelselectmetrics(:,1)));
    avgct = 0;
    for iii=1:framect-2,
      [ans,micsort] = sort(channelselectmetrics(:,iii));
      micsortbest = micsort(micct-locchannelct+1:micct);
      if ((max(squeeze(powers(ii,micsortbest,iii))) > powersratiothresh*maxpowers) & (min(squeeze(maxcorrs(ii,micsortbest,:,iii))) > maxcorrthresh) & (mindelaydiffs(ii,iii) < rejectlocthresh) & (coherence(ii,iii) > coherencethresh) ), %(mindelaylocs(ii,iii) > 0),
        if ((minlocs(ii,iii,1) >= lowcorners(whichsource,1)) & (minlocs(ii,iii,2) >= lowcorners(whichsource,2)) & (minlocs(ii,iii,1) <= highcorners(whichsource,1)) & (minlocs(ii,iii,2) <= highcorners(whichsource,2))),
          channelselectmean = channelselectmean + channelselectmetrics(:,iii);	% accumulate channel metrics for overall selection
          for iiii=(iii-1)*framesize-mindelay+1:fftblockstep:iii*framesize-mindelay-1,
            avgct = avgct + 1;
            %  FFT data, cross-correlate in each frequency bin, and sum into correlation estimate  (SEE FMV CODE)
            fftbuf = fft(xtemp(:,iiii:iiii+fftlen-1),[],2);
            for iiiii=fmindex(ii):fmaxindex(ii),
              Rxcorr(iiiii,:,:) = Rxcorr(iiiii,:,:) + reshape(fftbuf(:,iiiii)*fftbuf(:,iiiii)',1,micct,micct);
            end
          end
        end
      end
    end
avgct
    [junk,bestchannelsort] = sort(channelselectmean)
    refchannelindextrue = bestchannelsort(length(bestchannelsort))
    usechannels = sort(bestchannelsort(micct-extractchannelct+1:micct))
    for iii=1:extractchannelct,
      if (usechannels(iii) == refchannelindextrue), refchannelindex = iii; end
    end
    %  Extract steering vectors
    Hfsteer = zeros(extractchannelct,fftlen);
eigratio = zeros(1,fftlen);
    for iii=fmindex(ii):fmaxindex(ii),
      Rxcorruse = squeeze(Rxcorr(iii,usechannels,usechannels));
      [V,D] = eig(Rxcorruse);
      [eigvals,eigindex] = sort(diag(real(D)));
eigratio(iii) = eigvals(length(eigindex)-1)/eigvals(length(eigindex));
      Hfsteer(:,iii) = V(:,eigindex(length(eigindex)));
      Hfsteer(:,iii) = Hfsteer(:,iii)/Hfsteer(refchannelindex,iii);
    end
    %plot(real(ifft(Hfsteer.')));
    
    %  Recover desired source
    
    if (extractxyzflag(whichsource) == 1)	% Recover using the colocated directional array
      xtemp = filterit(micarray(startsamp:endsamp,5:8)',freqbands(ii,1),freqbands(ii,2),fs);
    else		% Use the spatially separated array to recover this frog
      xtemp = zeros(extractchannelct,size(xfilt,2));
      for iii=1:extractchannelct,
        xtemp(iii,1:xlen-meandelay(usechannels(iii))) = xfilt(usechannels(iii),meandelay(usechannels(iii))+1:xlen);
      end
    end
    
    %  apply FMV beamformer with steering vectors estimated for that source
    thissource(whichsource,:) = 2*real(robustafunc(xtemp',conj(Hfsteer),fftlen))';	% you only give it the positive-frequency steering vectors, so it comes out complex and half size
    

    
    figure(40+whichsource)
    plot(tottimes,real(thissource(whichsource,:)))
    grid on
    
  end
end


%   DONE!
