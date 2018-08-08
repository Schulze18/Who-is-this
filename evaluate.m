close all
clear V_z
clear Vk_z

N = length(signal);

zk = best_zk
Fk = best_Fk
Rz = best_Rz

N1 = best_N1
N2 = best_N2
% 
% 
% zk =...
% ...
%  [    1.880475972313699e-01...
%      9.365717710619779e-01...
%      9.847564547223975e-01...
% ]
% 
% Fk =[ ...
%      5.833780822381298e+00...
%      8.407584132419515e+00...
%      2.985734671496988e+02...
% ]
% 
% Rz =...
% ...
%      8.636353915427998e+00
%      
     
T=Ts;

for i=1:length(zk)
    Vk_z(i) = tf((1-2*abs(zk(i))*cos(2*pi*Fk(i)*T)+abs(zk(i))^2),[1 -2*abs(zk(i))*cos(2*pi*Fk(i)*T)  abs(zk(i))^2],Ts);
    
end

V_z = Rz*tf([1 -1 zeros(1,2*length(zk)-1)],1,Ts);


for i=1:length(zk)
    V_z = V_z*Vk_z(i);
end

entrada3 = glottal_model(entrada2,N1,N2);

saida_calc = lsim(V_z,entrada3);


%soundsc(saida_calc,Fs)
%soundsc(signal,Fs)

plot(saida_calc*max(signal)/max(saida_calc))
hold
plot(signal,'r')


xlabel('amostras')
ylabel('amplitude')
legend('sinal estimado', 'sinal original filtrado')




cd ..\Results


filename = 'ORIGINAL_FILTRADO_Bill_185-00003372.wav';
%filename = 'ORIGINAL_FILTRADO_Charlotte_3346-00006357.wav';
%filename = 'ORIGINAL_FILTRADO_Gale102-00005935.wav';
%filename = 'ORIGINAL_FILTRADO_Gale102-00005937.wav';
%filename = 'ORIGINAL_FILTRADO_Garth3092-00003917.wav';
%filename = 'ORIGINAL_FILTRADO_Garth3092-00007263.wav';
%filename = 'ORIGINAL_FILTRADO_Kelly1423-00004661.wav';
%filename = 'ORIGINAL_FILTRADO_Kelly1423-00007515.wav';
%filename = 'ORIGINAL_FILTRADO_Larry3143-00006878.wav';
%filename = 'ORIGINAL_FILTRADO_Larry3143-00006882.wav';
%filename = 'ORIGINAL_FILTRADO_Maria3017-00004490.wav';
%filename = 'ORIGINAL_FILTRADO_Moe_3125-00006208.wav';
%filename = 'ORIGINAL_FILTRADO_Moe_3125-00006209.wav';
audiowrite(filename,signal,Fs);



saida_calc = saida_calc*max(signal)/max(saida_calc);

filename = 'ESTIMADO_Bill_185-00003372.wav';
%filename = 'ESTIMADO_Charlotte_3346-00006357.wav';
%filename = 'ESTIMADO_Gale102-00005935.wav';
%filename = 'ESTIMADO_Gale102-00005937.wav';
%filename = 'ESTIMADO_Garth3092-00003917.wav';
%filename = 'ESTIMADO_Garth3092-00007263.wav';
%filename = 'ESTIMADO_Kelly1423-00004661.wav';
%filename = 'ESTIMADO_Kelly1423-00007515.wav';
%filename = 'ESTIMADO_Larry3143-00006878.wav';
%filename = 'ESTIMADO_Larry3143-00006882.wav';
%filename = 'ESTIMADO_Maria3017-00004490.wav';
%filename = 'ESTIMADO_Moe_3125-00006208.wav';
%filename = 'ESTIMADO_Moe_3125-00006209.wav';
audiowrite(filename,saida_calc,Fs);


soma_signals = (saida_calc(1:471000)/max(saida_calc)+signal(1:471000)'/max(signal))/2;


filename = 'COMPARACAO_Bill_185-00003372.wav';
%filename = 'COMPARACAO_Charlotte_3346-00006357.wav';
%filename = 'COMPARACAO_Gale102-00005935.wav';
%filename = 'COMPARACAO_Gale102-00005937.wav';
%filename = 'COMPARACAO_Garth3092-00003917.wav';
%filename = 'COMPARACAO_Garth3092-00007263.wav';
%filename = 'COMPARACAO_Kelly1423-00004661.wav';
%filename = 'COMPARACAO_Kelly1423-00007515.wav';
%filename = 'COMPARACAO_Larry3143-00006878.wav';
%filename = 'COMPARACAO_Larry3143-00006882.wav';
%filename = 'COMPARACAO_Maria3017-00004490.wav';
%filename = 'COMPARACAO_Moe_3125-00006208.wav';
%filename = 'COMPARACAO_Moe_3125-00006209.wav';
audiowrite(filename,soma_signals,Fs);

cd ..\07_12

