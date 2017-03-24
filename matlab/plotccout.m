%
%		plotccout.m
%		Douglas L. Jones
%		University of Illinois
%		July 6, 2007
%
%	plotcc.m: plots up to ten data vectors according to Ratnam's microphone color-code, with extensions for more channels:
%			Channel 1: cyan (yellow for Ratnam)
%			Channel 2: blue
%			Channel 3: green
%			Channel 4: red
%			Channel 5: black
%
%   NOTES:
%
%
function plotccout(data,ref,tottimes,fignum)

%   SETUP

figure(fignum)
channelct = min(size(data));

%
%  PLOT
%
hold off
plot(tottimes,ref,'k')
hold on
plot(tottimes,data(1,:),'c')
if (channelct > 1)
  plot(tottimes,data(2,:),'b')
  if (channelct > 2)
    plot(tottimes,data(3,:),'g')
    if (channelct > 3)
      plot(tottimes,data(4,:),'r')
      if (channelct > 4)
        plot(tottimes,data(5,:),'m')
        if (channelct > 5)
          plot(tottimes,data(6,:),'k')
          if (channelct > 6)
            plot(tottimes,data(7,:),'b-')
            if (channelct > 7)
              plot(tottimes,data(8,:),'k-')
              if (channelct > 8)
                plot(tottimes,data(9,:),'v')
                if (channelct > 9)
                  plot(tottimes,data(10,:),'p')
                end
              end
            end
          end
        end
      end 
    end    
  end
end
hold off
grid on

%  DONE

