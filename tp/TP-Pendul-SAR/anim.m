function [szOut]=anim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% anim.m: Fonction qui fait l'animation du pendule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% chargement des variables
load xtheta;
xtheta=ans';
xtheta= [xtheta.Time xtheta.Data];

% definition des limites de la figure d'animation
x1plus = max(xtheta(:,2));
x1moins = min(xtheta(:,2));
offset= 1;

% on dessine le premier
pose= [xtheta(1,2);xtheta(1,3)];
draw(pose,1);
axis([x1moins-offset x1plus+offset -0.6 0.6]);
refresh;

% on fait l'animation
cl0=clock;
for index=2:1:size(xtheta,1)
    while  etime(clock,cl0) < xtheta(index,1)
     etime(clock,cl0);
    end;
    pose= [xtheta(index,2);xtheta(index,3)];
    draw(pose,index);
 end;

function [szOut] = draw(pose,index)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% draw.m: Fonction qui dessine le pendule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global hc;
% Definition de la geometrie des corps et reperes
N= 32;
for i=1:1:N+1
  alpha= 2*(i-1)*pi/N;
  co(i)= cos(alpha);
  so(i)= sin(alpha);
end;

r_roue= 0.03; %rayon des roues
length= 0.4; %longueur de la barre
dx= 0.2; %demi-longueur du chariot
dy= 0.1; %hauteur du chariot

Roue= r_roue*[co(1) co(2) nan co(2) co(3) nan co(3) co(4) nan co(4) co(5) nan...
      co(5) co(6) nan co(6) co(7) nan co(7) co(8) nan co(8) co(9) nan...
      co(9) co(10) nan co(10) co(11) nan co(11) co(12) nan co(12) co(13) nan...
      co(13) co(14) nan co(14) co(15) nan co(15) co(16) nan co(16) co(17) nan...
      co(17) co(18) nan co(18) co(19) nan co(19) co(20) nan co(20) co(21) nan...
      co(21) co(22) nan co(22) co(23) nan co(23) co(24) nan co(24) co(25) nan...
      co(25) co(26) nan co(26) co(27) nan co(27) co(28) nan co(28) co(29) nan...
      co(29) co(30) nan co(30) co(31) nan co(31) co(32) nan co(32) co(33);
      so(1) so(2) nan so(2) so(3) nan so(3) so(4) nan so(4) so(5) nan...
      so(5) so(6) nan so(6) so(7) nan so(7) so(8) nan so(8) so(9) nan...
      so(9) so(10) nan so(10) so(11) nan so(11) so(12) nan so(12) so(13) nan...
      so(13) so(14) nan so(14) so(15) nan so(15) so(16) nan so(16) so(17) nan...
      so(17) so(18) nan so(18) so(19) nan so(19) so(20) nan so(20) so(21) nan...
      so(21) so(22) nan so(22) so(23) nan so(23) so(24) nan so(24) so(25) nan...
      so(25) so(26) nan so(26) so(27) nan so(27) so(28) nan so(28) so(29) nan...
      so(29) so(30) nan so(30) so(31) nan so(31) so(32) nan so(32) so(33)];
  
Chariot=  [-dx dx nan -dx dx nan -dx -dx nan dx dx; ...
    0 0 nan dy dy nan 0 dy nan 0 dy];
Barre= [0 0; 0 length];
Rep = [-100 100 nan 0 0; -2.2*r_roue -2.2*r_roue nan -0.02-2.2*r_roue 0.02-2.2*r_roue];

x= pose(1);
theta= pose(2);
orientation= [cos(-theta) -sin(-theta);sin(-theta) cos(-theta)];

for i=1:95
   Roue_Gauche(:,i)=  Roue(:,i)+ [x-0.75*dx;-r_roue];
   Roue_Droite(:,i)=  Roue(:,i)+ [x+0.75*dx;-r_roue];
end;   
for i=1:11
  Chariot(:,i)= Chariot(:,i)+[x;-r_roue];
end;
for i=1:2
  Barre(:,i)= [x;0]+ orientation*Barre(:,i);
end;

% draw the figure
if index==1   
% figure properties Initialization
  figure(1);
  clf;
  set(gcf, 'Name', 'Animation Window');
  set(gcf,'Position',[ 10 50 950 600]);
  set(1,'BackingStore','on');
  hold on;
  set(gca, 'UserData', hc,'NextPlot', 'add','Visible', 'on','DataAspectRatio', [1 1 1], ...
	  'Color', 'w', 'SortMethod', 'childorder');
  hc(1)=plot(Roue_Gauche(1,:),Roue_Gauche(2,:),'k-','Linewidth',3);
  hc(2)=plot(Roue_Droite(1,:),Roue_Droite(2,:),'k-','Linewidth',3);
  hc(3)=plot(Chariot(1,:),Chariot(2,:),'r-','Linewidth',3);
  hc(4)=plot(Barre(1,:),Barre(2,:),'b-','Linewidth',3);
  hc(5)=plot(Rep(1,:),Rep(2,:),'g-','Linewidth',3);
  drawnow;
  zoom(1.3);
  xlabel('x');
  ylabel('y');
  box on;
else
  set(hc(1),'XData',Roue_Gauche(1,:),'YData',Roue_Gauche(2,:));
  set(hc(2),'XData',Roue_Droite(1,:),'YData',Roue_Droite(2,:));
  set(hc(3),'XData',Chariot(1,:),'YData',Chariot(2,:));
  set(hc(4),'XData',Barre(1,:),'YData',Barre(2,:));
  set(hc(5),'XData',Rep(1,:),'YData',Rep(2,:));
  drawnow;
end;
