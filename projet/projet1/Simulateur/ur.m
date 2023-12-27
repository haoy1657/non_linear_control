%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fonction qui definit v1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [output]= ur(tps)
if tps <= 5
  freq= 1/2;
  vr= 6;
  wr=2*sin(2*pi*freq*(tps-15)+pi/2);
  elseif tps <= 10
  freq= 1/2;  
  vr= -5;
  wr=2*sin(2*pi*freq*(tps-15)+pi/2);
elseif tps <= 15
  vr= 5;
  wr=0.3;  
elseif tps <= 19
  vr= -4;
  wr=0.5;  
else
  vr= 1;
  wr= 0; 
end;

% % 
% if tps <= 5
%   freq= 1/2;
%   vr= 20;
%   wr=2*sin(2*pi*freq*(tps-15)+pi/2);
%   elseif tps <= 10
%   freq= 1/2;  
%   vr= -20;
%   wr=2*sin(2*pi*freq*(tps-15)+pi/2);
% elseif tps <= 15
%   vr= 20;
%   wr=1;  
% elseif tps <= 19
%   vr= -20;
%   wr=1.5;  
% else
%   vr= 1;
%   wr= 0; 
% end;


%output=[5;0];
output=[vr+2*rand;wr+0.1*rand];
