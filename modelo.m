close all
clear all
clc

ComPausa = 1;



if ComPausa
    
    load dados_com_pausa.mat
    N = length(pitch3);
    Ts = 1/Fs;

    %throat = iddata([signal(1:N)'],[glottal_out' amp3_filt'],Ts);
    throat = iddata([signal(1:N)'],[glottal_out'.*amp3_filt'],Ts);
else
    load dados_sem_pausa.mat
    Ts = 1/Fs;
    %throat = iddata([sinal_sem_pausa'],[pitch_sem_pausa' amp_sem_pausa'],Ts);
end

%throat.InputName  = {'Pitch';'Amplitude'};
%throat.InputName  = {'Pitch*Amplitude'};
%throat.OutputName = {'Voice'};

%mp = idss(throat(1:100000))

%m3 = armax(throat(100001:200000),'na',5,'nb',[5 5],'nc',2,'nk',[1 1]);
%m3 = armax(throat(100001:200000),'na',5,'nb',[5 ],'nc',5,'nk',[5 ]);

arx_opts = arxOptions('Focus','simulation');
oe_opts = oeOptions('Focus','simulation');
Gz_arx_ident = arx(throat(1:100000),[2,1,1],arx_opts);
Gz_oe_ident = oe(throat(1:100000),[2,1,1],oe_opts);

compare(throat, Gz_arx_ident,Gz_oe_ident)

%z = tf('z',Ts);


zk=.999;
Fk = 2;
T=1e-3;
Vk_z1 = tf((1-2*abs(zk)*cos(2*pi*Fk*T)+abs(zk)^2),[1 -2*abs(zk)*cos(2*pi*Fk*T)  abs(zk)^2],Ts);

zk=.999;
Fk = 10;
T=1e-3;
Vk_z2 = tf((1-2*abs(zk)*cos(2*pi*Fk*T)+abs(zk)^2),[1 -2*abs(zk)*cos(2*pi*Fk*T)  abs(zk)^2],Ts);


zk=.999;
Fk = 80;
T=1e-3;
Vk_z3 = tf((1-2*abs(zk)*cos(2*pi*Fk*T)+abs(zk)^2),[1 -2*abs(zk)*cos(2*pi*Fk*T)  abs(zk)^2],Ts);


Vk_z=Vk_z1*Vk_z2*Vk_z3;

saida_calc = lsim(Vk_z,glottal_out.*amp3_filt,(1:N)*Ts);

%sound(1/13*saida_calc,Fs)

%m4 = tfest(throat(1:end/2,'Voice',:),2,1,'Ts',Ts);

%compare(throat(1:100000),m3,m4,mp)

aux_a = max(saida_calc);
aux_b = max(signal);

mean((((saida_calc/aux_a)'-signal/aux_b))./((signal/aux_b)+1e-12))

%mean(((saida_calc).^2'-(signal.^2))./((signal+.00000000000001).^2))


%mean(((signal.^2)-((.9*signal).^2))./((signal+.00000000000001).^2))
mean(((signal)-((.9*signal)))./((signal+.00000000000001)))

figure
if 0
plot(saida_calc(100000:101000),'r')
hold
plot(signal(100000:101000))
else
plot(signal)
hold
plot(saida_calc,'r')
end
%compare(throat(1:end),m4)


%[YH, FIT, X0] = compare(throat(1:end),m4);

%sound(YH.OutputData,Fs)