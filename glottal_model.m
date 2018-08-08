function [ glottal_out ] = glottal_model( impulsetrain,N1,N2 )
factor=4;
n1 = 1:N1;
n2 = N1:N1+N2;
g = [.5*(1-cos(pi*n1/N1)) cos(pi*(n2-N1)/(2*N2))];

glottal_out = conv(g,impulsetrain);

end

