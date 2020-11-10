%%% Calcula e simula lei de controle GPC

function [y,u,Du,r] = gpc(B,A,nu,ny,Wu,Wy,ref,dist,noise)

sizey = size(A,1);
if size(B,2)==sizey;B=[B,zeros(sizey,sizey)];end

%%   Matrizes de predição

[H,P,Q] = predicao(A,B,ny);

%%   Lei de Controle

[Nk,Dk,Pr] = leimpc(H,P,Q,nu,Wu,Wy,sizey);

%Parametros de simulação
nNk = size(Nk,2)/sizey;
nDk = size(Dk,2)/sizey;
init = max([nNk,nDk])+2;
y = zeros(sizey,init);
u = y;
Du = u;
r = u;
d=u;
runtime = size(ref,2);

%% Simulaçao
for i=init:runtime-2;

%Atualiza lei de controle
d(1:sizey,i+1)=dist(:,i+1);
ypast = y(:, i:-1:i+1-nNk)+noise(:, i:-1:i+1-nNk);
Dupast = Du(:, i-1:-1:i-nDk) ;
upast = u(:, i-1);
rfut = ref(:,i+1); 

Dufast = Pr*rfut - Nk*ypast(:) - Dk*Dupast(:);
Du(:,i) = Dufast(1:sizey);
u(:,i)=u(:,i-1)+Du(:,i);

%Simula processo
upast2 = u(:,i:-1:i-nDk);
ypast2 = y(:, i:-1:i+2-nNk);
y(:,i+1) = -A(:,sizey+1:nNk*sizey)*ypast2(:) + B*[upast2(:)] + d(:,i+1);
r(:,i+1) = ref(:,i+1);

end

