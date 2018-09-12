figure
subplot(3,3,1) 
plot(hist(1,:),'*')
title(1)
ylim([0 0.5])
for i = 2:9
   subplot(3,3,i) 
   if 100*(i-1) > 750
        plot(hist(800,:),'*')
        title(800)
   else
        plot(hist(100*(i-1),:),'*')
        title(100*(i-1))
   end
   ylim([0 0.5])
   
end