y = rad2deg(y);
plot(t,y,'linewidth',2);
title('Time vs Theta1, Theta2, Theta1-dot, Theta2-dot')
lgd = legend('theta1','theta2','theta1-dot','theta2-dot');
lgd.FontSize = 14;

xlabel("Time in Seconds")
ylabel("Angle in Radians")

figure
subplot(2,2,1)      
plot(t,y(:,1),'linewidth',2);          
title('Time vs Theta 1')
xlabel("Time in Seconds")
ylabel("Angle in Radians")

subplot(2,2,2)       
plot(t,y(:,2),'linewidth',2);       
title('Time vs Theta 2')
xlabel("Time in Seconds")
ylabel("Angle in Radians")

subplot(2,2,3)       
plot(t,y(:,3),'linewidth',2);      
title('Time vs Theta 1-Dot')
xlabel("Time in Seconds")
ylabel("Angle in Radians")

subplot(2,2,4)       
plot(t,y(:,4),'linewidth',2);        
title('Time vs Theta 2-Dot')
xlabel("Time in Seconds")
ylabel("Angle in Radians")
