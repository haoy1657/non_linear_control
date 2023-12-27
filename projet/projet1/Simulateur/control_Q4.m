
function [output]=control_Q4(inp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fonction qui calcule la loi de commande
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% On renomme les variables
  pr= inp(1:2); % position de la r�f�rence 
  thetar= inp(3); % orientation de la r�f�rence 
  Xr=[pr;thetar]; % Etat complet de la r�f�rence
  p= inp(4:5); % position de la voiture
  theta= inp(6); % orientation de la voiture
  phi= inp(7); % angle de braquage de la voiture
  X= [p;theta;phi];
  vr= inp(8); % vitesse lin�aire de la r�f�rence
  wr= inp(9); % vitesse angulaire de la r�f�rence 
  
  % Distance entre les essieus avant et arri�re
  L=2;

  eps=0.3;

  % phir=0;
  phir=atan(L*vr*wr/(eps+vr^2));
% Calcul de la commande
  u= [0;0];
  % etude du premier sous-systeme, de premier ordre
  k1=1; % choix
  R_moinsthetar=[cos(Xr(3)) sin(Xr(3));-sin(Xr(3)) cos(Xr(3))];
  p_tilde=R_moinsthetar*(p-pr);
  theta_tilde=theta-thetar;
  phi_tilde=phi-phir;%on considere phir est nul,si il est pas nul ,on doit utiliser
  %une formule approximation pour le calculer 
  u(1)=-k1*p_tilde(1)+vr; %u_tilde(1)=u(1)-u_r1
% etude du deuxieme sous-systeme, de troisieme ordre
% on va essyer trouver les gain pertinent pour stabiliser ce systeme
  
  A=[0 vr 0;0 0 vr/L;0 0 0]; 
  B=[ 0 0 1]';
  pole_souhaite=[-2 -3 -5 ]; % choix
  K=place(A,B,pole_souhaite); % methode du placement de poles
  u_tilde2=-K*[p_tilde(2) theta_tilde phi_tilde]';
  u(2)=u_tilde2+wr;
  disp(k1);
  

 % On sort la commande, ainsi que l'�tat de la r�f�rence et de la voiture
 % pour visualisation gr�ce � la fonction anim
output= [u;X;Xr];

function out= rot(theta)
out= [cos(theta) -sin(theta); sin(theta) cos(theta)];
  