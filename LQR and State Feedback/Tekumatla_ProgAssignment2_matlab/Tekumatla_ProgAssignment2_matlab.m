clc;clear;close all;
%Programming Assignment 2 

%% Run this program first , once the gain values are calculated then run robot_control.m program 
% Manually updating gain values in the robot_control.m file is  not
% required
% All the variables or parameterized in the ODE function 
% If you wish to change the values such as mass , and lengths , just change
% them once at lines 41-49

%%

% a) Computing Equilibrium points
fprintf("------------Finding out  Equilibrium Points----------------\n");
syms l1 l2 r1 r2 theta1 theta2 I1 I2 m1 m2 theta1_dot theta2_dot g theta1_ddot theta2_ddot u1 u2 

EoM1 = theta1_ddot*(m2*l1^2 + 2*m2*cos(theta2)*l1*r2 + m1*r1^2 + m2*r2^2 + I1 + I2) - theta2_dot*(l1*m2*r2*theta1_dot*sin(theta2) + l1*m2*r2*sin(theta2)*(theta1_dot + theta2_dot)) - u1 + theta2_ddot*(m2*r2^2 + l1*m2*cos(theta2)*r2 + I2) - g*m2*r2*sin(theta1 + theta2) - g*l1*m2*sin(theta1) - g*m1*r1*sin(theta1);


EoM2 = theta2_ddot*(m2*r2^2 + I2) - u2 + theta1_ddot*(m2*r2^2 + l1*m2*cos(theta2)*r2 + I2) - g*m2*r2*sin(theta1 + theta2) + l1*m2*r2*theta1_dot*sin(theta2)*(theta1_dot + theta2_dot) - l1*m2*r2*theta1_dot*theta2_dot*sin(theta2);

EoM1 = subs(EoM1,[u1,u2,theta1_dot,theta2_dot,theta1_ddot,theta2_ddot],[0,0,0,0,0,0]); %Dynamics of the System are Zero

EoM2 = subs(EoM2,[u1,u2,theta1_dot,theta2_dot,theta1_ddot,theta2_ddot],[0,0,0,0,0,0]); %Dynamics of the System are Zero

eq_sol = solve(EoM1==0,EoM2==0,[theta1,theta2]);

disp("Equilibbrium Points Are");

disp("Theta1 , Theta2 ,Theta1_Dot, Theta2_Dot");
Equilibrium_Points = [eq_sol.theta1';eq_sol.theta2';[0,0,0];[0,0,0]]';
Equilibrium_Points(4,:) = [pi,pi,0,0]; %This is also a equilibrium point , But Matlab could not compute it . 

disp(Equilibrium_Points);

% disp("Theta2");
% disp(eq_sol.theta2);
%%
fprintf("------------Jacobian Linearization----------------\n");
% b) Jacobian Linearization
m1 = 1;
m2 = 1;
l1 = 1;
l2 = 1;
r1 = 0.45;
r2 = 0.45;
I1 = 0.084;
I2 = 0.084;
g  = 9.81;


X_dot = sym('X_dot' , [4,1]); %State Space Model 

X_dot(1) = theta1_dot;
X_dot(2) = theta2_dot;
X_dot(3) =  (I2*u1 - I2*u2 + m2*r2^2*u1 - m2*r2^2*u2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) + g*l1*m2^2*r2^2*sin(theta1) + I2*g*l1*m2*sin(theta1) + I2*g*m1*r1*sin(theta1) - l1*m2*r2*u2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*m1*m2*r1*r2^2*sin(theta1) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);
X_dot(4) = -(I2*u1 - I1*u2 - I2*u2 - l1^2*m2*u2 - m1*r1^2*u2 + m2*r2^2*u1 - m2*r2^2*u2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1^3*m2^2*r2*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) - g*l1^2*m2^2*r2*sin(theta1 + theta2) - I1*g*m2*r2*sin(theta1 + theta2) + g*l1*m2^2*r2^2*sin(theta1) + I2*g*l1*m2*sin(theta1) + I2*g*m1*r1*sin(theta1) + l1*m2*r2*u1*cos(theta2) - 2*l1*m2*r2*u2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + 2*l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) + l1^2*m2^2*r2^2*theta2_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + g*l1^2*m2^2*r2*cos(theta2)*sin(theta1) - g*m1*m2*r1^2*r2*sin(theta1 + theta2) + I1*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*m1*m2*r1*r2^2*sin(theta1) + 2*l1^2*m2^2*r2^2*theta1_dot*theta2_dot*cos(theta2)*sin(theta2) + l1*m1*m2*r1^2*r2*theta1_dot^2*sin(theta2) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2) + g*l1*m1*m2*r1*r2*cos(theta2)*sin(theta1))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);

A = jacobian(X_dot,[theta1,theta2,theta1_dot,theta2_dot]);

B = jacobian(X_dot, [u1,u2]);
disp("State Matrix After Jacobian Linearization");
disp(A);
disp("Input Matrix After Jacobian Linearization");
disp(B);
%%

% C) Stability Properties around each Equilibrium Point
fprintf("------------Stability Test Around All Equilibrium Points----------------\n");

for point = Equilibrium_Points'
    disp(point');
    u1 =0; u2 =0 ;
    A_eq = double(subs(A,[theta1,theta2,theta1_dot,theta2_dot],point'));
    B_eq = double(subs(B,[theta1,theta2,theta1_dot,theta2_dot],point'));
    disp("State Matrix at this Equilibrium Point is");
    disp(A_eq);
    disp("Input Matrix at this Equilibrium Point is");
    disp(B_eq);
    eig_vals = eig(A_eq);
    disp("Eigen Values are");
    disp(eig_vals);
    stable = all(real(eig_vals)<=0.0000001);
    if stable==1
        disp('The system is stable at the above equilibrium point')
    else
        disp('The system is NOT stable at the above equilibrium point')
    end 
    disp("_______________________________________________________________________________________________________________");
end

%%
% d) Controllability Around the Upward Configuration
fprintf("------------Controllability Test Around the upward Configurwation----------------\n");
disp("Upward Configurwation is [0,0,0,0]");

A_up = double(subs(A,[theta1,theta2,theta1_dot,theta2_dot],[0,0,0,0]));

B_up = double(subs(B,[theta1,theta2,theta1_dot,theta2_dot],[0,0,0,0]));

if rank(ctrb(A_eq,B_eq))==4
        disp("The system is  Controllable at the above equilibrium point")
else
        disp("The system is  NOT Controllable at the above equilibrium point")
end 
%%
% e) State-Feedback Control at upward equilibrium Point
fprintf("------------State Feedbback Design----------------\n");

lambda = [-1+2i,-1-2i,-2,-1.5];

Kn = place(A_up,B_up,lambda);

disp("Eigen Values placed at:")
disp(lambda);

disp("Calculated Gains at this point are:")
disp(Kn);
%%
% f)Feedback Control Law Implementaion 

fprintf("------------Feedback Control Law Implementation----------------\n");

disp("Plotting Results")

figure(1);

[t,y] = ode45(@(t,y) RR(t,y,Kn,m1,m2,l1,l2,r1,r2,I1,I2,g), [0,10],[pi/6,pi/4,0,0]);

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
ylabel("Angular Velocity in Radians per second ")

subplot(2,2,4)       
plot(t,y(:,4),'linewidth',2);        
title('Time vs Theta 2-Dot')
xlabel("Time in Seconds")
ylabel("Angular Velocity in Radians per second ")

%Calculate  forces/torques 

input1 = zeros(size(y(:,1)));
input2 = zeros(size(y(:,1)));


for i = 1:size(y,1)
    input1(i) = -Kn(1,:) *[y(i,1) ; y(i,2) ; y(i,3) ;y(i,4)];
    input2(i) = -Kn(2,:) *[y(i,1) ; y(i,2) ; y(i,3) ;y(i,4)];
    
end

%Plot Forces/Torques
figure
subplot(2,1,1)       
plot(t,input1,'linewidth',2);        
title('Time vs Joint -1 Torque')
xlabel("Time in Seconds")
ylabel("Torque(N-m) ")

subplot(2,1,2)       
plot(t,input2,'linewidth',2);        
title('Time vs Joint -2 Torque')
xlabel("Time in Seconds")
ylabel("Torque(N-m) ")


%%
%ODE Function 
function dz  = RR(t,z,Kn,m1,m2,l1,l2,r1,r2,I1,I2,g)

dz = zeros(4,1);
z = num2cell(z);

[theta1 , theta2 , theta1_dot , theta2_dot] = deal(z{:}) ; 

if abs(theta1) > 2*pi
    theta1 = mod(theta1, 2*pi);
end 

if abs(theta2) > 2*pi
    theta2 = mod(theta2, 2*pi);
    
end 

U =  -Kn *[theta1 ; theta2 ; theta1_dot ;theta2_dot];

u1 = U(1,:);
u2 = U(2,:);
 
dz(1) = theta1_dot; 
dz(2) = theta2_dot; 
dz(3) = (I2*u1 - I2*u2 + m2*r2^2*u1 - m2*r2^2*u2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) + g*l1*m2^2*r2^2*sin(theta1) + I2*g*l1*m2*sin(theta1) + I2*g*m1*r1*sin(theta1) - l1*m2*r2*u2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*m1*m2*r1*r2^2*sin(theta1) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);
dz(4) = -(I2*u1 - I1*u2 - I2*u2 - l1^2*m2*u2 - m1*r1^2*u2 + m2*r2^2*u1 - m2*r2^2*u2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1^3*m2^2*r2*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) - g*l1^2*m2^2*r2*sin(theta1 + theta2) - I1*g*m2*r2*sin(theta1 + theta2) + g*l1*m2^2*r2^2*sin(theta1) + I2*g*l1*m2*sin(theta1) + I2*g*m1*r1*sin(theta1) + l1*m2*r2*u1*cos(theta2) - 2*l1*m2*r2*u2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + 2*l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) + l1^2*m2^2*r2^2*theta2_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + g*l1^2*m2^2*r2*cos(theta2)*sin(theta1) - g*m1*m2*r1^2*r2*sin(theta1 + theta2) + I1*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*m1*m2*r1*r2^2*sin(theta1) + 2*l1^2*m2^2*r2^2*theta1_dot*theta2_dot*cos(theta2)*sin(theta2) + l1*m1*m2*r1^2*r2*theta1_dot^2*sin(theta2) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2) + g*l1*m1*m2*r1*r2*cos(theta2)*sin(theta1))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);
end

