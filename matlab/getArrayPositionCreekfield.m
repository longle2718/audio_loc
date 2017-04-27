function [array_pos, nMic, npair, pair_mic, c] = getArrayPositionCreekfield(...
    excode,isplot)

% array position
% excode: code of experiment
% 1: Creekfield day 1
% date: Jun4 2013
%
% 2: Creekfield day 2
% date: Jun5 2013
%

% Night of June 4, 2013 (Day 1)
% field orientation
%            |N x
%            |
%  W y       |          E
% water ----------- tree
%            |
%            |
%            |S
% Night of Jun 5, 2013 (Day 2)
% field orientation
%            | water y
%            |
%            |          x
% tree  ----------- tree
%            |
%            |
%            |tree

% -------------------------------------------------------------------------
if nargin < 2
    isplot = 0;
end
if nargin < 1
    array_pos_version = 1;
end

if excode == 1 %Jun4    
        % distance between 2 consecutive mics in x direction:
        d = 30 * 0.3048; % unit: meter
        % 1. use measurement data for array position
        % array_pos in feet
        array_pos = [ -130  -10  10;... %M1
            -100  0   10;... %M2
            -70   0   3 ;... %M3
            -40   -10 3 ;... %M4
            -10   -10 10;... %M5
            20    0   10;... %M6
            50    -10 3 ;... %M7
            80    0   3 ;... %M8
            ];
        % convert array_pos to meter
        array_pos = array_pos*0.3048;
        c = 347; 
elseif excode == 2 % Jun5
    % distance between 2 consecutive mics in x direction:
    d = 45 * 0.3048; % unit: meter
    % 1. use measurement data for array position % path is about 3' above
    % water level
    % % array_pos in feet
    array_pos = [ 225  9.5  10;... %M1
         180  -89/12   10 ;... %M2
         135  -13.5    3  ;... %M3
         90    0       2.5;... %M4
         45    0       10 ;... %M5        
         0     9       10;... %M6
        -45    0       3 ;... %M7
        -90    115/12  3 ;... %M8
        ];
    % angle between baseline 1 and 2
    theta = 2*asind((11.3/2)/45);
    % rotate baseline 2 clockwise an angle theta
    R = [cosd(theta) sind(theta); -sind(theta) cosd(theta)];
    array_pos(6:8,1:2) = array_pos(6:8,1:2)*(R');
    % convert array_pos to meter
    array_pos = array_pos*0.3048;    
end


nMic = size(array_pos,1);
npair = nchoosek(nMic,2);
pair_mic = zeros(npair,2);
count = 1;
for imic = 1:nMic-1
    for jmic = imic+1:nMic
        % name of mic in the pair
        pair_mic(count,1) = imic;
        pair_mic(count,2) = jmic;
        count = count + 1;
    end
end

% plot array configuration
if isplot == 1
    load('colorspec.mat');
    mycolor = colorspec(1:nMic,:);
    mysize = 80;
    mic_list = 1:nMic;
    mic_list_str = num2str(mic_list'); mic_text = cellstr(mic_list_str);
    dx1 = 0.1; dy1 = 0.1; dz1 = 0.1; % displacement from the plotted point
    figure(1);
    scatter3(array_pos(:,1), array_pos(:,2), array_pos(:,3),  mysize, mycolor,'filled');
    xlim([(min(array_pos(:,1)) - 5)  (max(array_pos(:,1)) + 5)]);
    ylim([(min(array_pos(:,2)) - 5)  (max(array_pos(:,2)) + 5)]);
    zlim([(min(array_pos(:,3)) - 1)  (max(array_pos(:,3)) + 1)]);
    text(array_pos(:,1) + dx1, array_pos(:,2) + dy1,array_pos(:,3) + dz1,mic_text);
    grid on;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('Array position');
    view(0,90);
    % axis equal;
end
                
            

