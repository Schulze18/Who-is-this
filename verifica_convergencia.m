z_fit = zeros(1,particle_number);
for i=1:particle_number
        
        zk = particles(1:2:(2*ordem_trato-1),i);
        Fk = particles(2:2:(2*ordem_trato),i);
        Rz = particles(end-2,i);
        N1 = round(particles(end-1,i));
        N2 = round(particles(end,i));
        [z, lag, Rz] = desired_function(zk, Fk, Rz, N1, N2, Ts, entrada, saida, MAX_LAG, goal);
       z_fit(i) = z;
end
figure(10)
hold on
plot(z_fit,'*')
ylim([0 1])
%%
z_fit = zeros(1,gwolfs_number)
for i=1:gwolfs_number
        
        zk = gwolfs(1:2:(2*ordem_trato-1),i);
        Fk = gwolfs(2:2:(2*ordem_trato),i);
        Rz = gwolfs(end-2,i);
        N1 = round(gwolfs(end-1,i));
        N2 = round(gwolfs(end,i));
        [z, lag, Rz] = desired_function(zk, Fk, Rz, N1, N2, Ts, entrada, saida, MAX_LAG, goal);
       z_fit(i) = z;
end

figure(10)
hold on
plot(z_fit,'*')
ylim([0 1])