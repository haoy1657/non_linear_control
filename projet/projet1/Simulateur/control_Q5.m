function [output]=control_Q5(inp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
  eps=0.1;
  L=2; % Distance entre les essieus avant et arri�re
  phir=atan(L*vr*wr/(eps+vr^2));
  R_moinsthetar=[cos(Xr(3)) sin(Xr(3));-sin(Xr(3)) cos(Xr(3))];
  p_tilde=R_moinsthetar*(p-pr);
  theta_tilde=theta-thetar;
  phi_tilde=phi-phir;
  
 
  % changement de variables
  xi1=p_tilde(1);
  xi2=p_tilde(2);
  xi3=tan(theta_tilde);
  xi4=tan(phi_tilde+phir)/(L*(cos(theta_tilde))^3);

% Calcul de la commande
  u= [0;0];
  k1=1; % choix
  w1=-k1*xi1+vr-wr*xi2; % loi de commande, voir rapport page 11
  u(1)=w1/cos(theta_tilde); % w1=u1*cos(theta_tilde)

  A=[0 1 0;0 0 1;0 0 0]; 
  B=[0 0 1]';
  pole_souhaite=[-4 -6 -8 ]; % choix
  %pole_souhaite=[-0.5+0.5i -0.5-0.5i -0.5 ];
  K=place(A,B,pole_souhaite);%k2 k3 k4


%xi2 xi2_dot xi2_dotdot
  xi2_dot=w1*xi3-wr*xi1;
  xi2_dotdot=w1^2*xi4-(w1*vr*(1+xi3^2)*tan(phir))/L;
%alpha(xi) et gamma(xi) -> voir rapport page 12
  alpha_xi=w1^2*((1+xi3^2)^3+L^2*xi4^2)/(L*(1+xi3^2)^1.5);
  gamma_xi=(3*w1^2*xi3*xi4)*(w1*xi4/(1+xi3^2)-vr*(tan(phir))/L)-(2*xi3*w1*vr*(tan(phir))/L)*(w1*xi4-vr*(1+xi3^2)*((tan(phir))/L));

%loi de commande w2, voir rapport page 12 
  w2=(-K(1)*xi2-K(2)*xi2_dot-K(3)*xi2_dotdot-gamma_xi)/alpha_xi;
  u(2)=w2;
 % On sort la commande, ainsi que l'�tat de la r�f�rence et de la voiture
 % pour visualisation gr�ce � la fonction anim
output= [u;X;Xr];

function out= rot(theta)
out= [cos(theta) -sin(theta); sin(theta) cos(theta)];
  