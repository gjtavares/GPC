close all; clc; clear all

%% Funcao de transferencia em Laplace

NumGs = 3.2132e-04;
DenGs = [3.4051 3.2132e-04];
Gs = tf(NumGs, DenGs);
stepinfo(Gs) % ts = 4.1457e+04

NumP1s = -0.0221;
DenP1s = [3.4051 3.2132e-04];
P1s = tf(NumP1s, DenP1s);

NumP2s = 2.3923e-07;
DenP2s = [3.4051 3.2132e-04];
P2s = tf(NumP2s, DenP2s);

%% Periodo de amostragem
Ts = 1382; % ts/15 <= Ts <= ts/6,  Ts = (ts/2)/15, ts/2 --> Tempo de acomodação em malha fechada

%% Discretizacao do modelo em tempo continuo:

Gz = c2d(Gs, Ts, 'zoh')
P1z = c2d(P1s, Ts, 'zoh')
P2z = c2d(P2s, Ts, 'zoh')



%% Modelo de G
DenGz=[1 -0.8777]; 
NumGz=[0 0.1223];
sizey=1;

%% Parametros
ny=10;  % Horizonte de predição
nu=5;   % Horizonte de Controle

per_amost = Ts;
t_final = Ts*120;
N = t_final/per_amost; %N° amostras
tempo = 0:per_amost:t_final-per_amost-1; %Vetor de tempo
t = 0:per_amost:t_final-3*per_amost; %Vetor de tempo

w = 35*ones(1,N);
wp =w(1,1:N-1);
%p=0.5*[zeros(1,30);ones(1,50);zeros(1,140)];
p=0.2*[zeros(1,40),ones(1,50),zeros(1,80)];
pp=p(1,1:N-1);

[y,u,Du,r] = gpc(NumGz,DenGz,nu,ny,1,1,w,p,[zeros(1,120)]);

%% Plot

figure(1)
subplot(2,1,1);
plot(tempo,wp)
hold on
plot(tempo,y)
grid on
legend('T1ref(t)','T1(t)')
title('Temperatura do Tanque')
xlabel('t (s)')
ylabel('Temperatura')
subplot(2,1,2);
hold on
plot(tempo,pp)
title('Pertubação')
xlabel('t(s)')
ylabel('Abertura da válvula')
grid on

figure(2)
plot(t,u)
title('Ação de Controle')
xlabel('t(s)')
ylabel('U')
grid on

