close all
clear all
%clc

% load handel
% 
% signal = filter(1/8*[1 1 1 1 1 1 1 1],1,y(1:70000)');
% 
% signal = y(1:70000)';

person = 2;
if person == 1
    load ../../signals/Originais/Bill_185-00003372.mat
    data2 = [data(10000:270371) data(283630:500000)]; % corta bips indesejados
    %no sinal Bill
elseif person ==2
    load ../../signals/Originais/Garth3092-00007263.mat
    data2 = data;
end

% soundsc(data2,44100)
%%

%load ../signals/Originais/Charlotte_3346-00006357.mat
%load ../signals/Originais/Gale102-00005935.mat
%data2 = data;

%Palavras

%sound(data2(1:24000),Fs)
%sound(data2(24000:52000),Fs)
%sound(data2(52000:59000),Fs)



%load ../signals/Originais/Charlotte_3346-00006357.mat
%load ../signals/Originais/Gale102-00005935.mat


janela_ms = 50;

M = 1000;%Fs*janela_ms*1e-3;
O = 4096;

OV = M/2; % window overlap

No = length(data2);

data2 = data2(1:floor(No/M)*M+1);

Fs = fs;

[B, A] = cheby1(3,.1,1/2*800/Fs);


signal = filter(B,A,data2);

%signal = filter(B,A,data(:,1)');


N = length(signal);



time2 = (1:M-OV:N-M)/Fs;

time3 = (1:N)/Fs;

%spectgm = zeros(N,M/2+1+O);
%periodograma = zeros(N,M+O);

correlograma = zeros(floor(N/(M-OV)),O);

%r = zeros((N-M)/(M-OV),M);
r = zeros(floor((N-M)/(M-OV)),2*M-1);

%h = spectrum.cov;
%h = spectrum.welch;
pitch2 = zeros(1,round((N-M-1)/(M-OV))+1);
amp2 = zeros(1,round((N-M-1)/(M-OV))+1);

amp3 = zeros(1,N);
pitch3 = zeros(1,N);
pitch4 = zeros(1,N);


tic
k=1;
for i = 1:M-OV:N-M
    %k = round((i-1)/(M-OV))+1;
    %r(k+1) = ([zeros(1,N+k) signal(i:i+M-1) zeros(1,N-k)]*[zeros(1,N-k) conj(signal(i:i+M-1)) zeros(1,N+k)]')/(N); %sequencia de autocovariancia
                                                                                        %standard biased ACS estimate
%      for L = 0:2*M-1                                                                             
%         r(k,L+1) = ([zeros(1,M+L) signal(i:i+M-1) zeros(1,M-L)]*[zeros(1,M-L) signal(i:i+M-1) zeros(1,M+L)]')/(M); %sequencia de autocovariancia
%      end                                                                                %standard biased ACS estimate
     
    media_pitch = mean(signal(i:i+M-1));
    r(k,:) = conv(signal(i:i+M-1)-media_pitch,signal(i+M-1:-1:i)-media_pitch)/(M);   %sequencia de autocovariancia
                                                            %standard biased ACS estimate
    
    correlograma(k,:) = abs(fft(r(k,:),O));
    %spectgm(k,:) = psd(h,[signal(i:i+M-1) zeros(1,O)],'Fs',Fs) ;
    %spectgm(k,:) = psd([signal(i:i+M-1) zeros(1,O)]) ;
    %[amp(k) pitch(k)] = max(spectgm(k,1:end/2));
     [amp2(k), pitch2(k)] = max(correlograma(k,1:end/2));
     
     amp2(k) = sqrt(amp2(k));
     pitch2(k)=(pitch2(k))*Fs/((M-OV)+O);
     
     if pitch2(k)<50 % elimina nível DC indesejado
         pitch2(k)=50;
     end
     
     %if k>1 && pitch2(k)==0,pitch2(k)=pitch2(k-1);end
         
     %pitch3(i:i+M-OV) = pitch3(i:i+M-OV) + (pitch2(k))/2;
     pitch3(i:i+M-OV) = (pitch3(i:i+M-OV) + pitch2(k))/2;
     pitch3(i+M-OV:i+M) = pitch2(k);
     
     

 
     pitch4(i) = pitch2(k);
     %amp4(i:i+M) = amp2(k);
     
     amp3(i:i+M-OV) = (amp3(i:i+M-OV) + amp2(k))/2;
     amp3(i+M-OV:i+M) = amp2(k);
     
     %sound(signal(i:i+M-1),Fs)
     k = k +1;
end
toc
mean(pitch2) 


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filtro no Sinal de Amplitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ordem_filtro = 20;
amp2_filt = filter(ones(1,Ordem_filtro)/Ordem_filtro,1,amp2);
amp3_filt = filter(ones(1,Ordem_filtro)/Ordem_filtro,1,amp3);

k=1;

nivel_amp_corte = mean(amp3_filt)/40;

% sinal_sem_pausa = signal(1);
% for index = 2:length(amp3_filt)
%     if amp3_filt(index)>= nivel_amp_corte
%         %corte(index)=1;
%         sinal_sem_pausa = [sinal_sem_pausa signal(index)];
%         %pitch_sem_pausa(k*(M-OV):k*(M-OV)+(M-OV)) = pitch3(index*(M-OV):index*(M-OV)+M-OV);
%         %amp_sem_pausa(k*(M-OV):k*(M-OV)+(M-OV)) = amp3_filt(index*(M-OV):index*(M-OV)+M-OV);
%         %k = k+1;
%     end
%     %der_corte(index) = corte(index)-corte(index);
%     
% end

% for index = 1:length(pitch_sem_pausa)
%     if index>1 && pitch_sem_pausa(index) == 0
%         %pitch_sem_pausa(index) = pitch_sem_pausa(index-1);
%     end
% end

%sinal_sem_pausa = filter(B,A,sinal_sem_pausa);

%sound(signal,Fs)

%subplot(311)
%plot(time2, pitch*Fs/(M+O)/2)
subplot(311)
plot(time3, pitch3,'.')
title('Sinais com pausa')
legend('pitch')
subplot(312)
plot(time3, pitch3/100,'g')
hold
plot((1:N)/Fs,signal)
 

legend('sinal','pitch')
subplot(313)
plot(time3, (amp3_filt)/10)
hold 
plot(time3,signal(1:length(time3)),'r')

legend('amplitude','sinal')



% figure
% 
% subplot(311)
% plot(pitch_sem_pausa,'.')
% title('Sinais sem pausa')
% legend('pitch')
% subplot(312)
% plot(sinal_sem_pausa)
% hold 
% plot(pitch_sem_pausa/10000,'g')
% legend('sinal','pitch')
% subplot(313)
% plot(amp_sem_pausa)
% hold 
% plot(sinal_sem_pausa,'r')
% 
% legend('amplitude','sinal')


pitch3 = filter(ones(1,(M-OV))/((M-OV)),1,pitch3);
for index=1:length(pitch3)
    if pitch3(index)<50
        pitch3(index)=50; % remove nivel DC indesejado
    end
end
figure
plot(time3, pitch3,'g')
legend('pitch')

impulsetrain = impulse_train_generator(pitch3*1,Fs);

N1 = 39*4;
N2 = 15*4;

glottal_out2 = glottal_model(impulsetrain,N1,N2);

glottal_out = glottal_out2(1:N);
hold on
plot(time3,impulsetrain,'r')

plot(time3,glottal_out,'y')


amp3_filt = filter(ones(1,(M-OV)/2)/((M-OV)/2)*1/15,1,amp3_filt);
figure
plot(time3, amp3_filt,'g')
legend('amplitude')



pitch_aux = pitch3;

pitch3 = pitch4;



%save('dados_sem_pausa.mat','pitch_sem_pausa', 'sinal_sem_pausa', 'amp_sem_pausa','Fs');

%save('dados_com_pausa.mat','pitch3', 'impulsetrain','glottal_out', 'signal', 'amp3_filt','Fs');
%sound(sinal_sem_pausa,Fs)

if person == 1
    save('dados_com_pausa_bill.mat','pitch3', 'impulsetrain','glottal_out', 'signal', 'amp3_filt','Fs');
elseif person == 2
    save('dados_com_pausa_garth.mat','pitch3', 'impulsetrain','glottal_out', 'signal', 'amp3_filt','Fs');
end
    
    