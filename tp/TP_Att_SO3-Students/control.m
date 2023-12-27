function [output]=control(inp)

% On d�finit les param�tres m�canique du corps
J= [1 -1 0; -1 2 1; 0 1 2];
B= [1 1 0.2;0.1 1 0; -0.1 -0.2 1];

% On renomme les variables
qs= inp(1);
qv= inp(2:4);
q=[qs; qv];
w= inp(5:7);

xi= 0.75;
puls= 1.5;
kq= 2*puls^2;
kw= 2*xi*puls;

% Solution Partie II-1
% control= zeros(3,1);
% ub = zeros(3,1);
% dqv = w/2;
% 
% for i=1:3
%     ub(i)=-2*xi*puls*dqv(i)-puls^2*qv(i);
% end
% 
% %Jb = 0.5 * J;
% control = 2 * inv(B) * J * ub;

% Jb=O.5*J-> SYST�ME INSTABLE (divergence), r�sultats obtenus
% insatisfaisants

%%Solution Partie II-2
%k= 2.3; % � peu pr�s proportionnel � puls^2
%K= 0.75 * eye(3,3); 
%control= (-inv(B) * k * qv * qs) - (2 * inv(B) * K * w);

%%Solution Partie II-3

K= 0.75 * eye(3); 
Kq = 2.5 * eye(3);
Kq=-Kq;

dqv= 0.5*(qs*w+cross(qv,w));

control= inv(B)*(-qs*qv+cross(w,J*w)+J*Kq*dqv-K*(w-Kq*qv));

output= [control;q;w];

end
