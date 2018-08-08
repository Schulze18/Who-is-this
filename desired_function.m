function [out, lag, Rz] = desired_function(zk, Fk, Rz, N1, N2, Ts, entrada, signal, MAX_LAG,goal)
format long e



%N = length(signal);

entrada = glottal_model(entrada,N1,N2);

%zk=.999;
%Fk = 2;
T=Ts;



%V_z = Rz*tf([1 -1 zeros(1,2*length(zk)-1)],1,Ts);
V_z = Rz*tf([1 -1 ],1,Ts);

Vk_z(1:length(zk)) = tf(0,1,Ts);

for i=1:length(zk)
    Vk_z(i) = tf((1-2*abs(zk(i))*cos(2*pi*Fk(i)*T)+abs(zk(i))^2),[1 -2*abs(zk(i))*cos(2*pi*Fk(i)*T)  abs(zk(i))^2],Ts);
    V_z = V_z*Vk_z(i);
end


if isstable(V_z)
    
    saida_calc = lsim(V_z,entrada);
else
    saida_calc = zeros(1,length(signal));
end

Rz = max(signal)/max(saida_calc);

%RMSerror = sqrt(mean((saida_calc'-signal).^2)/mean((signal).^2));
%RMSerror = sqrt(mean((saida_calc'-signal).^2));


% aux_X =xcorr(saida_calc/max(saida_calc),signal/max(signal),MAX_LAG);
% 
% [A, lag] = max(aux_X);
% 
% 
% 
% out = [1/A mean(abs(aux_X))];



aux_X = xcorr(saida_calc/max(saida_calc),signal/max(signal));
%aux_X =conv(saida_calc,signal);

A = max(aux_X);

lag=0;

%out = ([1/A]);
out = (goal*A);

%A = 100*(1-norm(signal-saida_calc')/norm(signal-mean(signal)));
%out = abs([1/A]);



%out = mean(abs(((saida_calc/aux_a)'./((signal/aux_b)+1e-15)-1)));
%out = mean(abs(((saida_calc/aux_a)'-signal/aux_b)));



