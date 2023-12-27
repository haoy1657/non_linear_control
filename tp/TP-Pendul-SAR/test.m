clear all;% Defintion des paramètres du système
m= 0.2; % masse de la barre
M= 1; % masse totale (chariot + barre)
l= 0.2; % demi-longueur de la barre
J=  m*l^2/3; % inertie de la barre
g= 9.81; % cste de gravité
mu= 0.005; % coeff. de frottement visqueux barre
tau= 0.2; % coeff. de frottement visqueux chariot

% Calcul des dérivées secondes
Delta_theta= M*(J+m*l^2)-m^2*l^2*(cos(0))^2;

A=              [0 0 1 0 
                 0 0 0 1 
                 0 -m^2*l^2*g/Delta_theta -tau*(J+m*l^2)/Delta_theta m*l*mu/Delta_theta
                 0 M*m*g*l/Delta_theta tau*m*l/Delta_theta -mu*M/Delta_theta];
B=[0 0 (J+m*l^2)/(M*(J+m*l^2)-m^2*l^2) -m*l/(M*(J+m*l^2)-m^2*l^2) ]';
disp(eig(A));
%il existe un pole dont la partie reel est positive donc , ce systeme est
%instable 
%methode lqr pour trouver le meilleur k 
Q=[100 0 0 0
    0 10 0 0
    0 0 1 0
    0 0 0 1 ]; 

R=  1;
%verifie commanbility 
cy=rank(ctrb(A,B));
disp('cy');
disp(cy);
K=lqr(A,B,Q,R);
%apres tester les bruit avec les frequence different, on constate que plus
%la frequence du bruit est grand ,moins l impact du bruit sur le systeme
%est petit
C=[1 0 0 0
    0 1 0 0];
D=0;
%verifier observabilite
oy=rank(obsv(A,C));
disp('oy');
disp(oy);
%为什么要添加观测器，状态空间方程的第一个方程是控制方程，我们默认其中所有的空间变量是可测的
%但是实际上方程的第二个方程告诉我们也就是矩阵c，只有几个变量是可观测的，比如我们这里
%可观变量是x和theta，[1 1 0 0],那么为了让第一个方程的控制可行，我们需要对不可观的变量
%进行预测，也就有了观测器
p=[-30+7.5i -30-7.5i  -40+12i -40-12i];
L=place(A',C',p);
L=L';

A=[0 1 0;0 0 1;0 0 0];
B=[0 ;0;1 ];
polesoute=[-1+2i -2-2i -5];
K=place(A,B,polesoute);
disp(K);
% u2_ksi=-K*[p2_ksi;theta_ksi;phi_ksi]
% u2=u2_ksi;
% u1_ksi=-1*p1_ksi;
% u1=u1_ksi+vr




