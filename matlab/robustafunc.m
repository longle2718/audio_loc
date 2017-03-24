%
%		robustafunc.m
%		Douglas L. Jones
%		University of Illinois
%		June 27, 2007
%
%	robustafunc.m: This code implements, as a function rather than a script,
%		the FMV beamformer with an
%		arbitrary number of sensors and arbitrary steering vectors
%		with a robustification that updates steering vectors adaptively within
%		a set maximum deviation from nominal to maximize output energy.
%		The idea is to track or learn the correct steering vectors
%
%	HISTORY
%	Modified from robusta.m
%
%	DEVELOPER NOTES
%  - Added compensation for non-unity window envelope in overlap-save some time ago
%	- Added tapered fade-in/fade-out window on 3/24/09, making it a mixed
%    overlap-save/add.  Did this because of ticking artifacts, which it cured
%
%   NOTES:
%       Should put damping factor in running recursive correlation calculation
%
%       Should zero-pad and process out past the end of the signal
%
%		Need to add anti-wraparound-artifact prevention, at least as option.
%		Any reason to taper it?  I can't think of one
%
%		Pass in the fmindex and fmaxindex, and save computations!
%
%
%function [xout] = robustafunc(x,evec,binct)
function [xout] = robustafunc(x,evec,binct)

%   SETUP

stopcircconvflag = 0;	%  Set to 1 to truncate filter to prevent any circular convolution wraparound
%binct = 1024;		%  number of FFT bins
binct2 = binct/2;
winlen = binct/2;	%  length of data window
corrlen = 2*binct;		%  number of samples to average in correlation estimate
corrlen2 = corrlen/2;
blklen = binct/8;		%  number of samples over which to apply same beamformer
blklen2 = blklen/2;

fadeblklen = 2*blklen;
fadeblklen2 = fadeblklen/2;
fadewindow = ones(1,fadeblklen);
fadewindow(1:fadeblklen2) = [0.5 - 0.5*cos(pi*[1:fadeblklen2]/fadeblklen2)];
fadewindow(fadeblklen-fadeblklen2+1:fadeblklen) = [0.5 + 0.5*cos(pi*[0:fadeblklen2-1]/fadeblklen2)];
fadewindow = fadewindow';


regfac = 0.05;		%  regularization factor
epsilon = 0.00000001;	%  divide-by-zero control

%   Create test signal

len = max(size(x));
channelct = min(size(x));
channelctsq = channelct*channelct;


%   Set parameters for robust adaptation

truesteerflag = 0;	% 1 means don't adapt steering vector, 0 adapts it (robust)
maxdev = 0.20;
maxdiff = 0.01;		% limit on magnitude of change vector per step relative to original magnitude

eveco = evec;

%   Initialize arrays

wb = hann(winlen);%hamming(winlen);
w = zeros(binct,channelct);	%  Make the STFT window
for ii=1:channelct,
  w(:,ii) = [zeros(1,(binct-winlen)/2) wb' zeros(1,(binct-winlen)/2)]'/max(wb);
%  w(:,ii) = hann(binct);
end
%w = ones(size(w));	% use a full-size boxcar window

wantiwrap = zeros(binct,channelct);
for ii=1:channelct,
  wantiwrap(:,ii) = fftshift([zeros(binct/4,1); ones(binct/2,1); zeros(binct/4,1)]);	%  Make the anti-circular-convolution window
end
%wantiwrap = ones(size(wantiwrap));


x = [ zeros(corrlen+blklen+binct,channelct); x; zeros(corrlen+blklen+binct,channelct) ];
xc = zeros(binct,channelctsq);

xout = zeros(max(size(x)),1);
xfout = zeros(binct,1);

fftbuf = zeros(binct,channelct);
corrbuf = zeros(binct,channelct^2);

%
%  APPLY ALGORITHM
%
weightsmat = zeros(binct,channelct);

epsilon = epsilon*max(max((abs(x))));

for iiii = 1:blklen:len,		% run through the data block by block
  iiii

  fftbuf = fft((w.*x(iiii+corrlen+blklen:iiii+corrlen+blklen+binct-1,:)));
  fftbufm = fft((w.*x(iiii+corrlen2+blklen:iiii+corrlen2+blklen+binct-1,:)));
  fftbufo = fft((w.*x(iiii:iiii+binct-1,:)));

%
%  Compute correlations and optimal filter in each frequency bin
%

  for iii=1:binct,
    xc(iii,:) = xc(iii,:) + reshape((fftbuf(iii,:)'*fftbuf(iii,:) - fftbufo(iii,:)'*fftbufo(iii,:)),1,channelctsq);	% recursively update running correlations in each bin
    
    Riii = reshape(xc(iii,:),channelct, channelct);
    Riii = Riii + (max(diag(Riii))*regfac + epsilon)*eye(channelct);           % normalized additive regularization of Riii

    Riiiinv = inv(Riii);
    if (truesteerflag == 1)	% standard Capon FMV
      newevec = evec(:,iii);
    else
% ROBUST EXTENSIONS FROM HERE

      thisevec = evec(:,iii);
      evecgrad = - 2*Riiiinv*thisevec/((thisevec'*Riiiinv*thisevec)^2 + eps);	% Compute gradient of output energy w.r.t. evec
      evecgradp = evecgrad - thisevec*(thisevec'*evecgrad)/(thisevec'*thisevec + eps);
      estep = - (evecgradp'*Riiiinv*thisevec + thisevec'*Riiiinv*evecgradp)/(2*evecgradp'*Riiiinv*evecgradp + epsilon);	% Compute step size to maximizing point
      if ( maxdiff*maxdiff*real(thisevec'*thisevec) < real(estep*estep*evecgradp'*evecgradp) )
        estep = sqrt( maxdiff*real(thisevec'*thisevec)/real(evecgradp'*evecgradp) );	% Limit to allowed range of deviation
      end
      newevec = thisevec + estep*evecgradp;	% Step to new steering vector
      diffvec = newevec - eveco(:,iii);
      if ( real(diffvec'*diffvec) > maxdev*maxdev*real(newevec'*newevec) )
        newevec = thisevec;
      end
      newevec = newevec*sqrt( (eveco(:,iii)'*eveco(:,iii))/(newevec'*newevec + eps) );	% Rescale to same magnitude as evec
      evec(:,iii) = newevec;
      
% TO HERE ...
    end

    weightsmat(iii,:) = (Riiiinv*newevec/(newevec'*Riiiinv*newevec + eps)).';         % compute optimal min-variance weights for the adjusted steering vector
%    xfout(iii) = fftbufm(iii,:)*weightsiii;            % apply filter in each frequency bin
  end

%
%   Compute time-domain output block
%
  if (stopcircconvflag == 1),
    weightsmattime = ifft(weightsmat,[],1).*wantiwrap;
    weightsmat = fft(weightsmattime,[],1);
  end
%whos
%pause

  xfout = sum(fftbufm.*weightsmat,2);
  xoutnow = ifft(xfout);
  xout(iiii+corrlen2+binct2+blklen-fadeblklen2:iiii+corrlen2+binct2+blklen+fadeblklen2-1) = xout(iiii+corrlen2+binct2+blklen-fadeblklen2:iiii+corrlen2+binct2+blklen+fadeblklen2-1) + fadewindow.*xoutnow(binct2-fadeblklen2+1:binct2+fadeblklen2)./w(binct2-fadeblklen2+1:binct2+fadeblklen2,1);

end

xout = xout(corrlen+blklen+binct+1:corrlen+blklen+binct+len);

%  DONE

