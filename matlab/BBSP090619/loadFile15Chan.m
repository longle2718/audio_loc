% PXI 4472, 8 chan
% Buffer size 500000 samples per chan
% Analog read and write to disk at 0.5*BuffSize = 250000 samples per chan
% File size = 1000000 samples per chan
% Data format: I32
% Range: -10V to +10V

fname = 'HcChorus_090619_001_16C_0143.dat';

PRECISION = 'int32';
minv = -10;
maxv = 10;

bytes_persample = 4; % double. Specify 'float64' in fread
nsamples_perblock = 250000;
nchans = 15;
dirInfo = dir(fname);
fsize = dirInfo.bytes;
nblocks = fsize/(nsamples_perblock*bytes_persample*nchans);
micarray = zeros(nblocks*nsamples_perblock, nchans);

fid = fopen(fname, 'r', 'native');

for j = 1:nblocks
    for i = 1:nchans
        [x, count] = fread(fid, nsamples_perblock, PRECISION);
        if (count ~= nsamples_perblock)
            error('Mysterious error: number of samples or blocks do not match')
        else
            fprintf('Block: %4d, Chan: %2d, samples read: %d\n', j, i, count);
            micarray((((j-1)*nsamples_perblock)+1):(j*nsamples_perblock),i) = x;
        end
    end
end
clear x
fclose(fid);

% scale
%minr = intmin(PRECISION);
%maxr = intmax(PRECISION);
minr = -2147483648;
maxr = 2147483647;



micarray = uint32todouble(micarray, minr, maxr, minv, maxv);
