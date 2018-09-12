clc, clear all, close all%, format long E
tic

%Define qual pessoa esta sendo analisada
person = 1;

woman=0;
multiplier=1;
if woman
    multiplier = 1.3
end

if person == 1
    load dados_com_pausa_bill.mat
    N = length(pitch3);
    Ts = 1/Fs;
elseif person == 2
    load dados_com_pausa_garth.mat
    N = length(pitch3);
    Ts = 1/Fs;
end

%%
%BILL

%%Trecho 1
% ini = 4000;
% fim = 16000;

%Trecho 2
% ini = 67000;
% fim = 76000;

%Trecho 3
% ini = 431000;
% fim = 440000;

% %Trecho 4
% ini = 4000;
% fim = 16000;

%Trecho 5
ini = 428000;
fim = 440000;

%Garth

% Trecho 1
% ini = 94000;
% fim = 106000;

% Trecho 2
% ini = 436000;
% fim = 448000;

%Charolote

% ini = 5000;
% fim = 13000;

%entrada = amp3_filt(ini:fim).*glottal_out(ini:fim);

entrada = amp3_filt(ini:fim).*impulsetrain(ini:fim);
saida = signal(ini:fim);


%entrada2 = amp3_filt.*glottal_out;
entrada2 = amp3_filt.*impulsetrain;

%Parametros da Simulação
N = 2; %Número de interações
gwolfs_number = 2000;
ordem_trato = 2;

NVz = 2*ordem_trato;

%field x vs y
%bottom = [0.7, 0, 0.7, 0, 0.7, 0]; % zk e Fk
%top = [1-1e-15, 20, 1-1e-15, 200, 1-1e-15, 10];


%bottom = zeros(1,NVz+1); % zk e Fk e Rz

%Inicialização dos limites superiores e inferiores de cada parametro
bottom = [];
top = [];
lim_F = 0;
for index=1:2:NVz
    %bottom = [bottom 0.7, lim_F];
    %bottom = [bottom .99+.001*index, 0];
    
    bottom = [bottom 0.7, lim_F];
    top = [top 1-1e-15, (1+index)*200];
    
    lim_F = top(end);
    
    %bottom = [bottom 0.7, 0];
    %top = [top 1-1e-15, 2000];
    
    %top = [top 1-1e-15, 2000];
    %top = [top 1-1e-15, index];
end
%%
% top = [1-1e-15 800 1-1e-15 1800]; %Teste apenas
% bottom = [0.70 0 0.70 800];       %Teste 

top = [top 40 40*2/multiplier 20*2/multiplier]; %Rz N1 N2
bottom = [bottom 1 4 2];

MAX_LAG = 1000;

goal = [(1/max(xcorr(saida/max(saida),saida/max(saida))))] % goal absoluto
eps = 0.47;
error_parameters = length(goal);
%%

%Parameteros SWO
% inertia = 0.9;
% initial_inertia = .9;
% self_confidence = 0.8;
% toward_best = 0.01;
% speed = .2;

%Parametros GWO
z_alfa = 0;
z_beta = 0;
z_delta = 0;
X_alfa = ones(1,length(bottom));
X_beta = ones(1,length(bottom));
X_delta = ones(1,length(bottom));
a = 2;
A = 1;
C = 1;

%Inicialização GWO
gwolfs_parameters = length(bottom);
for k=1:gwolfs_parameters
    gwolfs(k,:) = bottom(k)*ones(1,gwolfs_number)+(top(k)-bottom(k))*rand(1,gwolfs_number); 
end
%%
%Teste fixando os ganhos
% gwolfs(5,:) = 3;
% gwolfs(6,:) = 80;
% gwolfs(7,:) = 40;
% top((length(top)-2):length(top)) = [40 44*2/multiplier 22*2/multiplier]; %Rz N1 N2

%%
% for i = 1:gwolfs_number 
%      
%      zk = gwolfs(1:2:end-3,i);
%      Fk = gwolfs(2:2:end-3,i);
%      Rz = gwolfs(end-2,i);
%      N1 = round(gwolfs(end-1,i));
%      N2 = round(gwolfs(end,i));

     %Determina z,lag e Rz
%      [z, lag, Rz] = desired_function(zk, Fk, Rz, N1, N2, Ts, entrada, saida, MAX_LAG, goal);
%      
%      if z > z_delta
%          if z > z_beta
%              if z > z_alfa
%                  X_delta = X_beta;
%                  z_delta = z_beta;
%                  X_beta = X_alfa;
%                  z_beta = z_alfa;
%                  X_alfa = [z; zk; Rz; N1; N2; lag];
%                  z_alfa = z;
%              else
%                  X_delta = X_beta;
%                  z_delta = z_beta;
%                  X_beta = [z; zk; Rz; N1; N2; lag];
%                  z_beta = z;
%              end
%          else
%              X_delta = [z; zk; Rz; N1; N2; lag];
%              z_delta = z;
%          end
%      end
%end




%best_gwolf = 1;
%best_z = eps/100000000000000;
%z_best = ones(gwolfs_number,error_parameters)/100000000000000;
%gwolfs_best = 10*ones(gwolfs_parameters, gwolfs_number);



vetor_alfa = zeros(1,N);
vetor_beta = zeros(1,N);
vetor_delta = zeros(1,N);
for m=1:N
    hist(:,m) = gwolfs(:,1);
    m;
    %Verifica valor de fitness e atualiza Alfa, Beta e Delta
    for i=1:gwolfs_number
        
%         zk = gwolfs(1:2:end-3,i);
%         Fk = gwolfs(2:2:end-3,i);
%         Rz = gwolfs(end-2,i);
%         N1 = round(gwolfs(end-1,i));
%         N2 = round(gwolfs(end,i));
        
        zk = gwolfs(1:2:(2*ordem_trato-1),i);
        Fk = gwolfs(2:2:(2*ordem_trato),i);
        Rz = gwolfs(end-2,i);
        N1 = round(gwolfs(end-1,i));
        N2 = round(gwolfs(end,i));
        
        

        %Determina z,lag e Rz do Wolf atual
        [z, lag, Rz] = desired_function(zk, Fk, Rz, N1, N2, Ts, entrada, saida, MAX_LAG, goal);
        
        %Verifica verifica se é o novo Alfa, Beta ou Delta
        if z > z_delta(1)
            if z > z_beta(1)
                if z > z_alfa(1)
                    X_delta = X_beta;
                    z_delta = z_beta;
                    X_beta = X_alfa;
                    z_beta = z_alfa;
                    X_alfa = gwolfs(:,i);
                    z_alfa = [z; zk; Fk; Rz; N1; N2; lag];
                else
                    X_delta = X_beta;
                    z_delta = z_beta;
                    X_beta = gwolfs(:,i);
                    z_beta = [z; zk; Fk; Rz; N1; N2; lag];
                end
            else
                 X_delta = gwolfs(:,i);
                 z_delta = [z; zk; Fk; Rz; N1; N2; lag];
            end
        end
    end    

    %Armazena Alfa, Beta e Delta
    vetor_alfa(m) = z_alfa(1);  
    vetor_beta(m) = z_beta(1);
    vetor_delta(m) = z_delta(1);
    
    %Atualiza posição dos gwolfs
    for i=1:gwolfs_number
        %Distancia Alfa
        r1 = rand(length(X_alfa),1);r2 = rand(length(X_alfa),1);
        A1 = 2*a*r1-a; C1 = 2*r2;
        D_alfa = abs(C1.*X_alfa-gwolfs(:,i));
        X1 = X_alfa - A1.*D_alfa;

        %Distancia Beta
        r1 = rand(length(X_alfa),1);r2 = rand(length(X_alfa),1);
        A2 = 2*a*r1-a;
        C2 = 2*r2;
        D_beta = abs(C2.*X_beta-gwolfs(:,i));
        X2 = X_beta - A2.*D_beta;

        %Distancia Delta
        r1 = rand(length(X_alfa),1);r2 = rand(length(X_alfa),1);
        A3 = 2*a*r1-a;
        C3 = 2*r2;
        D_delta = abs(C3.*X_delta-gwolfs(:,i));
        X3 = X_delta - A2.*D_delta;

        gwolfs(:,i) = (X1+X2+X3)/3; 
        
        %Correção caso saia os limites definidos
        for j = 1:gwolfs_parameters
            while gwolfs(j,i) > top(j) || gwolfs(j,i) < bottom(j)
                if gwolfs(j,i) > top(j)
                  gwolfs(j,i) = top(j) - (gwolfs(j,i)-top(j));
                elseif gwolfs(j,i) < bottom(j)
                  gwolfs(j,i) = bottom(j) - (gwolfs(j,i)-bottom(j));
                end
            end
        end
    end
    %Fixa R,N1,N2
    %gwolfs(5,:) = 3; gwolfs(6,:) = 80; gwolfs(7,:) = 40;
    
    a = a - 2/N;
   
end

    z_alfa(1)
    z_beta(1)
    z_delta(1)
%%

%Calculo e Plot da resposta utilizando o melhor resultado - Alfa
m=0;
error = [1e3];
bla=0;
toc

%movie(movie1,1,2)
best_error = z_alfa(1);
%[zk, Fk] = particles(:,best_particle)

% %z_alfa = [z; zk; Fk; Rz; N1; N2; lag];
% zk = z_alfa(2);
% Fk = z_alfa(3);
% Rz = z_alfa(4);
% N1 = z_alfa(5);
% N2 = z_alfa(6);

zk = z_alfa(2:1:(2+ordem_trato-1));
Fk = z_alfa((2+ordem_trato):1:(2+2*ordem_trato-1));
Rz = z_alfa(end-3);
N1 = z_alfa(end-2);
N2 = z_alfa(end-1);


T=Ts;
%V_z = Rz*tf([1 -1 zeros(1,2*length(zk)-1)],1,Ts);
V_z = Rz*tf([1 -1 ],1,Ts);

for i=1:length(zk)
    Vk_z(i) = tf((1-2*abs(zk(i))*cos(2*pi*Fk(i)*T)+abs(zk(i))^2),[1 -2*abs(zk(i))*cos(2*pi*Fk(i)*T)  abs(zk(i))^2],Ts);
    V_z = V_z*Vk_z(i);
end


%%
entrada2 = glottal_model(entrada2,N1,N2);

saida_calc = lsim(V_z,entrada2);


%soundsc(saida_calc,Fs)
%soundsc(signal,Fs)
%%
plot(saida_calc/max(saida_calc))
hold
plot(signal/max(signal),'r')

legend('saida calc.', 'original')
%%

figure,
plot(xcorr(saida_calc(ini:fim)/max(saida_calc(ini:fim)),signal(ini:fim)/max(signal(ini:fim))))
hold
plot(xcorr(signal(ini:fim)/max(signal(ini:fim)),signal(ini:fim)/max(signal(ini:fim))))
legend('saida calc.', 'original')

figure,
plot(xcorr(saida_calc/max(saida_calc),signal/max(signal)))
hold
plot(xcorr(signal/max(signal),signal/max(signal)))
legend('saida calc.', 'original')

saida_calc = lsim(V_z,entrada);


%soundsc(saida_calc,Fs)
%soundsc(signal,Fs)

figure
plot(saida_calc/max(saida_calc))
hold
plot(saida/max(saida),'r')

legend('saida calc.', 'original')


toc
