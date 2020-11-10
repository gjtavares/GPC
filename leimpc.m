%%% Lei de Controle
%%  

function  [Nk,Dk,Pr,S,X,Prlong] =leimpc(H,P,Q,nu,Wu,Wy,sizey)


%Horizonte de controle
P1 = H(:,1:sizey*nu);

%% Matrizes de ponderação

WY=Wy;
WU=Wu;
L = eye(sizey);
npred = size(P1,1)/sizey;
for i = 2:npred;
   v=(i-1)*sizey+1:i*sizey; 
   WY(v,v) = Wy;
   WU(v,v) = Wu;
   L = [L;eye(sizey)];
 end
 WU = WU(1:nu*sizey,1:nu*sizey);
 

S = P1'*WY*P1 + WU; S=(S+S')/2;
X = [P1'*WY*P,P1'*WY*Q,-P1'*WY];
   

%% Parametros Lei de controle
M = inv(S);
Nk = M*P1'*WY*Q;
Dk = M*P1'*WY*P;
Prlong = M*P1'*WY;

Pr = Prlong*L;
X = [P1'*WY*P,P1'*WY*Q,-P1'*WY*L];

