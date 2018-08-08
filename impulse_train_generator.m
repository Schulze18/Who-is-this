function [ impulse_train ] = impulse_train_generator( pitch, Fs )


N = length(pitch);
impulse_train = zeros(1,N);

i=1;
while i<N
    
    impulse_train(i) = 1;
    i = i + round(Fs/pitch(i));

end

