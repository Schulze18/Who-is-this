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

%BILL

%Trecho 1
% ini = 4000;
% fim = 16000;

%Trecho 2
ini = 428000;
fim = 440000;



%Charolote

% ini = 5000;
% fim = 13000;

%entrada = amp3_filt(ini:fim).*glottal_out(ini:fim);

entrada = amp3_filt(ini:fim).*impulsetrain(ini:fim);
saida = signal(ini:fim);


%entrada2 = amp3_filt.*glottal_out;
entrada2 = amp3_filt.*impulsetrain;


%Parametros da Simulação
particle_number = 2000;
N = 10; %Número de interações
ordem_trato = 2;

NVz = 2*ordem_trato;

%field x vs y
%bottom = [0, 0, 0, 0, 0, 0]; % zk e Fk
%top = [1-1e-12, 20, 1-1e-12, 200, 1-1e-12, 10];


%bottom = zeros(1,NVz+1); % zk e Fk e Rz

bottom = [];
top = [];
lim_F = 0;
for index=1:2:NVz
   
    bottom = [bottom 0.7, lim_F];
    
    %top = [top 1-1e-15, 10*8^(index/2-0.5)];
    top = [top 1-1e-15, (1+index)*200];
    lim_F = top(end);
    
    
end

top = [top 40 40*2/multiplier 20*2/multiplier]; %Rz N1 N2
bottom = [bottom 1 4 2];

MAX_LAG = 1000;

goal = [(1/max(xcorr(saida/max(saida),saida/max(saida))))] % goal absoluto
eps = 0.70;
error_parameters = length(goal);

%Parameteros SWO
inertia = 5;
initial_inertia = 5;
self_confidence = 10;
toward_best = 10;
speed = 0.1;
z_best_global = 0;
X_best_global = zeros(length(bottom),1);
z_best_individual = zeros(particle_number,1);
X_best_individual = zeros(length(bottom)+2,particle_number);
particle_parameters = length(bottom);
index_best_global = 1;


%Inicialização do SWARM
for k=1:particle_parameters
    particles(k,:) = bottom(k)*ones(1,particle_number)+(top(k)-bottom(k))*rand(1, particle_number); 
    
end

for k = 1:particle_number
    sign_rand = (rand > 0.5)*2 - 1;
    particles_dpos(:,k) = sign_rand*initial_inertia*(bottom' + (top-bottom)'.*rand(particle_parameters,1));
end

% best_particle = 1;
% best_z = eps/100000000000000;
% z_best = ones(particle_number,error_parameters)/100000000000000;
% particles_best = 10*ones(particle_parameters, particle_number);


% m=0;
% error = [1e3];
% 
% bla=0;
% toc

%%
m = 0;
vetor_best_global = zeros(1,N);
%while max(abs(error) > eps)
%while ((z_best_global < eps) && (m < N) )
while m <N
    m = m+1
    %error_acum = 0;
    

    hist(:,m) = particles(:,1);
        
    for i=1:particle_number
        
        zk = particles(1:2:(2*ordem_trato-1),i);
        Fk = particles(2:2:(2*ordem_trato),i);
        Rz = particles(end-2,i);
        N1 = round(particles(end-1,i));
        N2 = round(particles(end,i));

        %best_z = desired_function(best_Q, best_R);
        [z, lag, Rz] = desired_function(zk, Fk, Rz, N1, N2, Ts, entrada, saida, MAX_LAG, goal);
        
        particles(end-2,i) = Rz; % atualiza o Rz da partícula recalculado dentro da desired_function
        
        if z > z_best_individual(i)
            z_best_individual(i) = z;
            X_best_individual(:,i) = [z; zk; Fk; Rz; N1; N2; lag];
            
            if z > z_best_global
                z_best_global = z;
                X_best_global = [z; zk; Fk; Rz; N1; N2; lag];
                index_best_global = i;
            end
        end

    end
    
    for i=1:particle_number
        %if particles(:,i) ~= particles(:,best_particle)
        if i ~= index_best_global
            for j=1:particle_parameters

                particles_dpos(j,i) = speed*(toward_best*rand*(X_best_global(j+1) - particles(j,i)) + self_confidence*rand*(X_best_individual(j+1,i) - particles(j,i)) + inertia*particles_dpos(j,i));
                particles(j,i) = particles_dpos(j,i) + particles(j,i);
                
                while particles(j,i) > top(j)  || particles(j,i) < bottom(j)
                    %bla = bla+1
                    if particles(j,i) > top(j)
                        %particles(j,i) = top(j)-abs(particles_dpos(j,i)/2);
                        particles(j,i) = top(j)-(particles(j,i)-top(j));
                        particles_dpos(j,i) = -particles_dpos(j,i);       
                    elseif particles(j,i) < bottom(j)
                        particles(j,i) = bottom(j)+(bottom(j)-particles(j,i));
                        particles_dpos(j,i) = -particles_dpos(j,i);
                    end
                end
            end
        end
    end
    vetor_best_global(m) = X_best_global(1);
end

%%

z_best_global
zk = X_best_global(2:1:(2+ordem_trato-1));
Fk = X_best_global((2+ordem_trato):1:(2+2*ordem_trato-1));
Rz = X_best_global(end-3);
N1 = X_best_global(end-2);
N2 = X_best_global(end-1);

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

figure
plot(vetor_best_global)

toc