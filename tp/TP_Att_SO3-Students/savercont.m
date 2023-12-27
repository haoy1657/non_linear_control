function [cont]=control(inp,methode,gainhomo,gainpt,gaintf1,gaintf2)

% Variables

g0=inp(1); g1=inp(2); g2=inp(3); g3=inp(4); g=[g0;g1;g2;g3];
v1=inp(5); v2=inp(6); v3=inp(7); v=[v1;v2;v3];
gr0=inp(15); gr1=inp(16); gr2=inp(17); gr3=inp(18); gr=[gr0;gr1;gr2;gr3];
vr1=inp(19); vr2=inp(20); vr3=inp(21); vr=[vr1;vr2;vr3];
dvr1=inp(22); dvr2=inp(23); dvr3=inp(24); dvr=[dvr1;dvr2;dvr3];
d2vr1=inp(25); d2vr2=inp(26); d2vr3=inp(27); d2vr=[d2vr1;d2vr2;d2vr3];
alf0= inp(28); tps=alf0; dalf0=1;
alf= inp(29);

% definitions generiques

  a=-1;
  S=[0 -1;1 0];
  gt= prodg(ginv(gr),g);
  Xgt= xg(gt);
  Adgtm= ad(ginv(gt));
  vt=v-Adgtm*vr;
  dgAdgtmvr= dgadgmv(gt,vr);
  dgt= Xgt*vt;
  dAdgtmvr= dgAdgtmvr*dgt+Adgtm*dvr;

if methode == 1  %%%%%%  Methode 1: Position tracking

  % gains de commande 
  
  L=[1; 0];
  Minv=inv([[1;0] S*L]);
  Kp= [gainpt(1) 0; 0 gainpt(2)]; 
  Kv= [gainpt(3) 0; 0 gainpt(4)];
  
  pb= [gt(1); gt(2)] + rot(gt(3))*L;
  dpb= rot(gt(3))*([vt(1);vt(2)]+vt(3)*S*L);
  temp= Minv*(-vt(3)*S*([vt(1);vt(2)]+vt(3)*S*L)-[v2*v3;-v1*v3]+...
          [dAdgtmvr(1);dAdgtmvr(2)]+dAdgtmvr(3)*S*L+...
          rot(-gt(3))*(-Kp*pb-Kv*dpb));
  u1= temp(1);
  u2= temp(2);
  dalf=0; 
  cont= [u1;u2;g;v;gt;gr;vr;dalf0;dalf];  

elseif methode == 2  %%%%%%  Methode 2: Fixed point Homogene  (ms97)

  % gains de commande 
  k1= gainhomo(1); k2= gainhomo(2); k3= gainhomo(3); k4= gainhomo(4);
  k5= gainhomo(5); omega= gainhomo(6);

  rho=(gt(2)^4+gt(3)^4+gt(4)^2+v3^2)^(1/4);
  u1=  -k1*(v1+k3*gt(2)-sin(omega*tps)*(k5*gt(4)+2*k5*v3)/rho);
  u2= -k2*(v2+k4*gt(3)+sign(a)*sin(omega*tps)*rho);
  dalf=0;
  cont= [u1;u2;g;v1;v2;v3;gt(2);gt(3);gt(4);gr;vr;dalf0;dalf];  

elseif methode == 3  %%%%%%  Methode 3: T.F. dh1=X1(h1) 

  eps1= gaintf1(1); eps2= gaintf1(2); k= gaintf1(3);
  k1= gaintf1(4); k2= gaintf1(5); k3= gaintf1(6);

  f1= eps1*sin(alf); f2= eps2*cos(alf); f3= a*eps1*eps2/4*sin(2*alf);
  daf1= eps1/eps2*f2; daf2= -eps2/eps1*f1; daf3= a*eps1*eps2/2*cos(2*alf);
  Minv= 2/a/eps1/eps2*[daf3 daf2; a*f1 1];
  daM= [0 f2; -a*daf1 -4*f3];
  h1= [cos(f1/2);sin(f1/2);0;0]; h1m= ginv(h1);
  Adh1= ad(h1); Adh1m= ad(h1m);
  Tf1= [1 0 0;0 1 0; 0 a*f1 1]; 
  
  gb= prodg(gt,h1m); gbm= ginv(gb);
  Xgb= xg(gb); 
  vstar= [-k1*gb(2)/gb(1);-k2*gb(3)/gb(1);-k3*gb(4)/gb(1)];
  dgvstar= [k1*gb(2)/gb(1)^2, -k1/gb(1), 0, 0;...
            k2*gb(3)/gb(1)^2, 0, -k2/gb(1), 0;...
            k3*gb(4)/gb(1)^2, 0, 0, -k3/gb(1)];
  Adgbm= ad(gbm);
  Advstar= Adgbm*vr+vstar; 
  dgAdvstar= dgadgmv(gb,vr)+dgvstar;
  xib= [v2-f2;v3+f3-a*f1*v2]-[0 1 0;0 0 1]*Advstar;
  R23= Xgb*Adh1*(Tf1*([0;xib]+Advstar)+[0;f2;f3]-Adgtm*vr);
  Z1= Xgb*Adh1*[1;0;0];
  temp= Minv*(-k*xib+[0;-a*v2*Advstar(1)]+...
              [0 1 0;0 0 1]*(Adgbm*dvr+dgAdvstar*R23));
  u2= temp(1);
  dalf= temp(2);
  df1= daf1*dalf; df2= daf2*dalf; df3= daf3*dalf;
  dM=dalf*daM;
  xib1= v1-df1-Advstar(1);
  dgb= xib1*Z1+R23;
  dAdvstar= dgAdvstar*dgb+Adgbm*dvr;
  dAdgbmdvr= Adgbm*d2vr+dgadgmv(gb,dvr)*dgb;
  dxib= [u2-df2;a*v1*v2+df3-a*df1*v2-a*f1*u2]-[0 1 0;0 0 1]*dAdvstar;
  dXgb= dxg(gb,dgb); 
  dAdh1= dad1(f1,df1);
  dR23= (dXgb+Xgb*dAdh1*Adh1m)*Adh1*(Tf1*([0;xib]+Advstar)+...
          [0;f2;f3]-Adgtm*vr)+...
       Xgb*Adh1*([0;0;a*df1*(xib(1)+Advstar(2))]+...
                  Tf1*([0;dxib]+dAdvstar)+[0;df2;df3]-dAdgtmvr);
  ddgAdvstar=  ddgadvstar(gb,vr,dgb,dvr);
  d2alf= [0 1]*Minv*(-dM*temp-k*dxib+...
             [0;-a*u2*Advstar(1)-a*v2*dAdvstar(1)]+...
             [0 1 0;0 0 1]*(dAdgbmdvr+ddgAdvstar*R23+dgAdvstar*dR23));
  d2f1= daf1*d2alf-f1*dalf^2;
  u1= dAdvstar(1)+d2f1-k*xib1+...
         xib'*([0; -a*v2]+[0 1 0;0 0 1]*dgAdvstar*Z1);
cont= [u1;u2;g;xib(1);xib(2);xib1;gt(2);gt(3);gt(4);gr;vr;dalf0;dalf];  

elseif methode == 4  %%%%%%  Methode 4: T.F. dh1= X2(h1)

  eps1= gaintf2(1); eps2= gaintf2(2); k= gaintf2(3);
  k1= gaintf2(4); k2= gaintf2(5); k3= gaintf2(6);

  f1= eps1*sin(alf); f2= eps2*cos(alf); f3= a*eps1*eps2/4*sin(2*alf);
  daf1= eps1/eps2*f2; daf2= -eps2/eps1*f1; daf3= a*eps1*eps2/2*cos(2*alf);
  Minv= 2/a/eps1/eps2*[daf3 daf2; a*f1 1];
  daM= [0 f2; -a*daf1 -4*f3];
  h1= [cos(f1/2);0;sin(f1/2);0]; h1m= ginv(h1);
  Adh1= ad(h1); Adh1m= ad(h1m);
  Tf1= [1 0 0;0 1 0; a*f1 0 1]; 
  
  gb= prodg(gt,h1m); gbm= ginv(gb);
  Xgb= xg(gb); 
  vstar= [-k1*gb(2)/gb(1);-k2*gb(3)/gb(1);-k3*gb(4)/gb(1)];
  dgvstar= [k1*gb(2)/gb(1)^2, -k1/gb(1), 0, 0;...
            k2*gb(3)/gb(1)^2, 0, -k2/gb(1), 0;...
            k3*gb(4)/gb(1)^2, 0, 0, -k3/gb(1)];
  Adgbm= ad(gbm);
  Advstar= Adgbm*vr+vstar; 
  dgAdvstar= dgadgmv(gb,vr)+dgvstar;
  xib= [v1-f2;v3+f3-a*f1*v1]-[1 0 0;0 0 1]*Advstar;
  R13= Xgb*Adh1*(Tf1*([xib(1);0;xib(2)]+Advstar)+[f2;0;f3]-Adgtm*vr);
  Z2= Xgb*Adh1*[0;1;0];
  temp= Minv*(-k*xib+[0;-a*v1*Advstar(2)]+...
              [1 0 0;0 0 1]*(Adgbm*dvr+dgAdvstar*R13));
  u1= temp(1);
  dalf= temp(2);
  df1= daf1*dalf; df2= daf2*dalf; df3= daf3*dalf;
  dM=dalf*daM;
  xib2= v2-df1-Advstar(2);
  dgb= xib2*Z2+R13;
  dAdvstar= dgAdvstar*dgb+Adgbm*dvr;
  dAdgbmdvr= Adgbm*d2vr+dgadgmv(gb,dvr)*dgb;
  dxib= [u1-df2;a*v1*v2+df3-a*df1*v1-a*f1*u1]-[1 0 0;0 0 1]*dAdvstar;
  dXgb= dxg(gb,dgb); 
  dAdh1= dad2(f1,df1);
  dR13= (dXgb+Xgb*dAdh1*Adh1m)*Adh1*(Tf1*([xib(1);0;xib(2)]+Advstar)+...
          [f2;0;f3]-Adgtm*vr)+...
       Xgb*Adh1*([0;0;a*df1*(xib(1)+Advstar(1))]+...
                 Tf1*([dxib(1);0;dxib(2)]+dAdvstar)+[df2;0;df3]-dAdgtmvr);
  ddgAdvstar=  ddgadvstar(gb,vr,dgb,dvr);
  d2alf= [0 1]*Minv*(-dM*temp-k*dxib+...
             [0;-a*u1*Advstar(2)-a*v1*dAdvstar(2)]+...
             [1 0 0;0 0 1]*(dAdgbmdvr+ddgAdvstar*R13+dgAdvstar*dR13));
  d2f1= daf1*d2alf-f1*dalf^2;
  u2= dAdvstar(2)+d2f1-k*xib2+...
         xib'*([0; -a*v1]+[1 0 0;0 0 1]*dgAdvstar*Z2);
cont= [u1;u2;g;xib(1);xib(2);xib2;gt(2);gt(3);gt(4);gr;vr;dalf0;dalf];  

end;


function [gp]=prodg(g1,g2) % group product of g1 and g2
  
  g1p= [g1(2);g1(3);g1(4)];
  g2p= [g2(2);g2(3);g2(4)];
  gp= [g1(1)*g2(1)-g1p'*g2p;g1(1)*g2p+g2(1)*g1p+vect(g1p,g2p)];

function [gi]= ginv(g) % group inverse of g

  gi= [g(1);-g(2);-g(3);-g(4)];

function [vf]= xg(g) % matrix X(g)  

  vf= 1/2*[-g(2) -g(3) -g(4); g(1) -g(4) g(3); g(4) g(1) -g(2);...
           -g(3) g(2) g(1)];

function [dvf]= dxg(g,dg) % value of d/dt X(g)

  dvf= 0.5*dg(1)*[0 0 0;1 0 0; 0 1 0; 0 0 1]+...
       0.5*dg(2)*[-1 0 0; 0 0 0; 0 0 -1; 0 1 0]+...
       0.5*dg(3)*[0 -1 0; 0 0 1; 0 0 0; -1 0 0]+...
       0.5*dg(4)*[0 0 -1; 0 -1 0; 1 0 0; 0 0 0];

function [adj]= ad(g) % Adjoint mapping

 g0=g(1); g1=g(2); g2=g(3); g3=g(4);
 adj= eye(3)+2*[-g2^2-g3^2, g1*g2-g0*g3, g1*g3+g0*g2;...
                g1*g2+g0*g3, -g1^2-g3^2, g2*g3-g0*g1;...
                g1*g3-g0*g2, g2*g3+g0*g1, -g1^2-g2^2];

function [dadj]= dad1(f,df) % time differential of Ad(h1) for dh1=X1(h1)

  dadj= df*[0 0 0; 0 -sin(f) -cos(f); 0 cos(f) -sin(f)];

function [dadj]= dad2(f,df) % time differential of Ad(h1) for dh1=X2(h1)

  dadj= df*[-sin(f) 0 cos(f); 0 0 0; -cos(f) 0 -sin(f)];

function [mat]=dgadgmv(g,v) % differential/g of Ad(g^{-1})v
  
  g0=g(1); g1=g(2); g2=g(3); g3=g(4); v1=v(1); v2=v(2); v3=v(3);
  mat= 2*[g3*v2-g2*v3, g2*v2+g3*v3, -2*g2*v1+g1*v2-g0*v3,...
          -2*g3*v1+g0*v2+g1*v3;...
          g1*v3-g3*v1, g2*v1-2*g1*v2+g0*v3, g1*v1+g3*v3,...
          -g0*v1-2*g3*v2+g2*v3;...
          g2*v1-g1*v2, g3*v1-g0*v2-2*g1*v3, g0*v1+g3*v2-2*g2*v3,...
          g1*v1+g2*v2];
          

function [mat]=ddgadvstar(g,v,dg,dv) 
     % time differential of the differential/g of Ad(g^{-1})v

  v1=v(1); v2=v(2); v3=v(3);
  mat =dgadgmv(g,dv)+...
       2*dg(1)*[0 0 -v3 v2; 0 v3 0 -v1; 0 -v2 v1 0]+...
       2*dg(2)*[0 0 v2 v3; v3 -2*v2 v1 0; -v2 -2*v3 0 v1]+...
       2*dg(3)*[-v3 v2 -2*v1 0; 0 v1 0 v3; v1 0 -2*v3 v2]+...
       2*dg(4)*[v2 v3 0 -2*v1; -v1 0 v3 -2*v2; 0 v1 v2 0];

function [ro]=vect(v,v2) % vector product

  ro= [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0]*v2; 
