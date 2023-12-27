function [szOut]=anim

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% anim.m: Fonction qui fait l'animation du corps rigide en rotation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% On charge les variables

load corps;
corps= ans';
corps= [corps.Time corps.Data];
load control;
control= ans';
control= [control.Time control.Data];

% on dessine le premier

rotr=rot([1;0;0;0]);
rots=rot([corps(1,2);corps(1,3);corps(1,4);corps(1,5)]);
co=[control(1,2);control(1,3)];
draw(rotr,rots,co,1);
axis(1.5*[-1.5 1.5 -1.5 1.5 -1.5 1.5]);
refresh;

% on fait l'animation
cl0=clock;
for index=2:1:size(corps,1)
    while  etime(clock,cl0) < corps(index,1)
     etime(clock,cl0);
    end;
    rotr=rot([1;0;0;0]);
    rots=rot([corps(index,2);corps(index,3);corps(index,4);corps(index,5)]);
    co=[control(index,2);control(index,3)];
    draw(rotr,rots,co,index);
 end;
 disp('Taper sur une touche pour fermer la fenetre');
 pause;
 close;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [szOut] = draw(rotr,rots,co,index)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% draw.m: Fonction qui dessine le repère de référence et le corps rigide 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global hc;

% Définition de la géométrie des corps et repères
% Corps contrôlé (forme conique)

N= 32;
for i=1:1:N+1
  alpha= 2*(i-1)*pi/N;
  co(i)= cos(alpha);
  so(i)= sin(alpha);
end;

rc= 1; %extrémités de la base du cone
lc= 1; %demi-longueur du cone

topc= rc/2*[co(1) co(2) nan co(2) co(3) nan co(3) co(4) nan co(4) co(5) nan...
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
      so(29) so(30) nan so(30) so(31) nan so(31) so(32) nan so(32) so(33);
      2*lc 2*lc nan 2*lc 2*lc nan 2*lc 2*lc nan 2*lc 2*lc nan...
      2*lc 2*lc nan 2*lc 2*lc nan 2*lc 2*lc nan 2*lc 2*lc nan...
      2*lc 2*lc nan 2*lc 2*lc nan 2*lc 2*lc nan 2*lc 2*lc nan...
      2*lc 2*lc nan 2*lc 2*lc nan 2*lc 2*lc nan 2*lc 2*lc nan...
      2*lc 2*lc nan 2*lc 2*lc nan 2*lc 2*lc nan 2*lc 2*lc nan...
      2*lc 2*lc nan 2*lc 2*lc nan 2*lc 2*lc nan 2*lc 2*lc nan...
      2*lc 2*lc nan 2*lc 2*lc nan 2*lc 2*lc nan 2*lc 2*lc nan...
      2*lc 2*lc nan 2*lc 2*lc nan 2*lc 2*lc nan 2*lc 2*lc];
  
  downc= rc*[co(1) co(2) nan co(2) co(3) nan co(3) co(4) nan co(4) co(5) nan...
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
      so(29) so(30) nan so(30) so(31) nan so(31) so(32) nan so(32) so(33);
      -lc -lc nan -lc -lc nan -lc -lc nan -lc -lc nan...
      -lc -lc nan -lc -lc nan -lc -lc nan -lc -lc nan...
      -lc -lc nan -lc -lc nan -lc -lc nan -lc -lc nan...
      -lc -lc nan -lc -lc nan -lc -lc nan -lc -lc nan...
      -lc -lc nan -lc -lc nan -lc -lc nan -lc -lc nan...
      -lc -lc nan -lc -lc nan -lc -lc nan -lc -lc nan...
      -lc -lc nan -lc -lc nan -lc -lc nan -lc -lc nan...
      -lc -lc nan -lc -lc nan -lc -lc nan -lc -lc];
  
  sidec= [-rc/2 -rc nan rc/2 rc nan 0 0 nan 0 0;
           0 0 nan 0 0 nan -rc/2 -rc nan rc/2 rc;
           lc -lc nan lc -lc nan lc -lc nan lc -lc];
       
% Repère lié au corps contrôlé
XY  = 1/1*[ 0 1.5 nan 1.5 1.5 nan 0 0 nan 0 0 nan 0 0 nan 0 0 ...
        nan 0 0 nan -0 0 nan 0 0;...
       0 0 nan 0 0 nan 0 0 nan 0 1.5 nan 1.5 1.5 nan 1.5 0 ...
        nan 0 0 nan 0 0 nan 0 0;...
        0 0 nan 0 0 nan 0 0 nan 0 0 nan -0 0 nan 0 0 ...
        nan 0 1.5 nan 1.5 1.5 nan 0 0];  
    
% Repère de référence
XYr = [ 0 1.5 nan 1.5 1.5 nan 0 0 nan 0 0 nan 0 0 nan 0 0 ...
        nan 0 0 nan -0 0 nan 0 0;...
       0 0 nan -0 0 nan 0 0 nan 0 1.5 nan 1.5 1.5 nan 1.5 0 ...
        nan 0 0 nan 0 0 nan 0 0;...
        0 0 nan 0 0 nan 0 0 nan 0 0 nan -0 0 nan 0 0 ...
        nan 0 1.5 nan 1.5 1.5 nan 0 0];
for i=1:26
  XY(:,i) = rots * XY(:,i);
end;
for i=1:95
  topc(:,i) = rots * topc(:,i);
  downc(:,i) = rots * downc(:,i);
end;
for i=1:11
  sidec(:,i) = rots *sidec(:,i);
end;
for i=1:26
  XYr(:,i) = rotr * XYr(:,i);
end;

% On dessine les objets
if index==1   
% Initialisation des propriétés de la figure
  figure(1);
  clf;
  set(gcf, 'Name', 'Animation Window');
  set(gcf,'Position',[ 10 50 950 600]);
  set(1,'BackingStore','on');
  hold on;
  set(gca, 'UserData', hc,'NextPlot', 'add','Visible', 'on','DataAspectRatio', [1 1 1], ...
	  'Color', 'w', 'SortMethod', 'childorder');
  view(3);
  hc(1)=plot3(XY(1,:),XY(2,:),XY(3,:),'b-','Linewidth',1);
  hc(2)=plot3(topc(1,:),topc(2,:),topc(3,:),'k-','Linewidth',3);
  hc(3)=plot3(downc(1,:),downc(2,:),downc(3,:),'g-','Linewidth',3);
  hc(4)=plot3(sidec(1,:),sidec(2,:),sidec(3,:),'k-','Linewidth',3);
  hc(5)=plot3(XYr(1,:),XYr(2,:),XYr(3,:),'r-','Linewidth',1);
  drawnow;
  box on;
  grid on;
else
  set(hc(1),'XData',XY(1,:),'YData',XY(2,:),'ZData',XY(3,:));
  set(hc(2),'XData',topc(1,:),'YData',topc(2,:),'ZData',topc(3,:));
  set(hc(3),'XData',downc(1,:),'YData',downc(2,:),'ZData',downc(3,:));
  set(hc(4),'XData',sidec(1,:),'YData',sidec(2,:),'ZData',sidec(3,:));
  set(hc(5),'XData',XYr(1,:),'YData',XYr(2,:),'ZData',XYr(3,:));
  drawnow;
  box on;
  grid on;
end;


function [mat]=rot(g)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcule la rotation associée à un quaternion  g
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 g0=g(1); g1=g(2); g2=g(3); g3=g(4);
 mat= eye(3)+2*[-g2^2-g3^2, g1*g2-g0*g3, g1*g3+g0*g2;...
                g1*g2+g0*g3, -g1^2-g3^2, g2*g3-g0*g1;...
                g1*g3-g0*g2, g2*g3+g0*g1, -g1^2-g2^2];
