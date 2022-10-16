clc;clear;close all;
syms theta1 theta2 theta1_dot theta2_dot a0 a1 a2 a3 b0 b1 b2 b3 t

%% 3.2(a) - Cubic Polynomial Tragectory 
q1 = a0+ a1*t + a2*t^2 + a3*t^3; % Cubic Polynomial for 1st joint

q2 = b0+ b1*t + b2*t^2 + b3*t^3; % Cubic Polynomial for 2nd joint

q_dot = jacobian([q1,q2],t);     % Time Derivative 

t0=0; tf =10; %Initial and final times 

eq1 = subs(q1,[t],[t0]); %a0+ a1*t0 + a2*t0^2 + a3*t0^3
eq2 = subs(q1,[t],[tf]); %a0+ a1*tf + a2*tf^2 + a3*tf^3
eq3 = subs(q_dot(1),[t],[t0]); 
eq4 = subs(q_dot(1),[t],[tf]);
eq5 = subs(q2,[t],[t0]); %b0+ b1*t0 + b2*t0^2 + b3*t0^3
eq6 = subs(q2,[t],[tf]); %b0+ b1*tf + b2*tf^2 + b3*tf^3
eq7 = subs(q_dot(2),[t],[t0]);
eq8 = subs(q_dot(2),[t],[tf]);

sol = solve([eq1==pi,eq2==0,eq3==0,eq4==0,eq5==pi/2,eq6==0,eq7==0,eq8==0],[a0,a1,a2,a3,b0,b1,b2,b3]);

cubic_1 = subs(q1,[a0,a1,a2,a3],[sol.a0,sol.a1,sol.a2,sol.a3]);
cubic_2 = subs(q2,[b0,b1,b2,b3],[sol.b0,sol.b1,sol.b2,sol.b3]);
disp("Cubic Polynomial for 1st Joint")
disp(cubic_1)
disp("Cubic Polynomial for 2nd Joint")
disp(cubic_2)
qd = [cubic_1;cubic_2];
qd_dot = jacobian(qd,t);
qd_ddot = jacobian(qd_dot,t);

%% 3.2(b) Manipulator Equation Form 
m1 = 1; m2= 1; l1=1;l2=1;r1=0.45;r2=0.45; I1 = 0.084; I2 = 0.084; g =9.81;
syms theta1 theta2
M = [m2*l1^2+2*m2*cos(theta2)*l2*r2+m1*r1^2+m2*r2^2+I1+I2   m2*r2^2+m2*r2*l1*cos(theta2)+I2;
    m2*r2^2+m2*r2*l1*cos(theta2)+I2                          m2*r2^2+I2];

C = [0                                         l1*m2*r2*theta1_dot*sin(theta2)+l1*m2*r2*sin(theta2)*(theta1_dot+theta2_dot);
    l1*m2*r2*(theta1_dot+theta2_dot)*sin(theta2)-l1*m2*r2*theta2_dot*sin(theta2)            0];
G = [-g*l1*m2*sin(theta1)-g*m1*r1*sin(theta1)-m2*g*r2*sin(theta1+theta2);
    -g*m2*r2*sin(theta1+theta2)];
%% Feedback Linearization 
%Virtual Control law  v = -[kp kd][q;q_dot];
A_new = [0 0 1 0;
         0 0 0 1;
         0 0 0 0;
         0 0 0 0];
B_new = [0 0;
         0  0;
         1 0;
         0 1];

lambda = [-3-2i,-3+2i,-3,-2.5];
disp("Eigen Values for virtual control input are placed at");
disp(lambda)
K = place(A_new,B_new,lambda);
disp("Gains are")
disp(K)
fprintf("------------Computing the ODE Outputs----------------\n");
[t,y] = ode45(@(t,y) RR(t,y,K,M,C,G,m1,m2,l1,l2,r1,r2, I1 , I2 ,g,qd,qd_dot,qd_ddot), [0,10],[deg2rad(200),deg2rad(125),0,0]);
theta1_desired =zeros(length(t),1);
theta2_desired =zeros(length(t),1);
theta1_dot_desired =zeros(length(t),1);
theta2_dot_desired =zeros(length(t),1);
u1                 =zeros(length(t),1);
u2                 =zeros(length(t),1);

for i=1:length(t)
    qd_sub = subs(qd,t(i));
    qd_dot_sub = subs(qd_dot,t(i));
    qd_ddot_sub = subs(qd_ddot,t(i));
    theta1_desired(i) = qd_sub(1);
    theta2_desired(i) = qd_sub(2);
    theta1_dot_desired(i) = qd_dot_sub(1);
    theta2_dot_desired(i) = qd_dot_sub(2);  
    M_n = double(subs(M,[theta1,theta2,theta1_dot,theta2_dot],[y(i,1), y(i,2),y(i,3),y(i,4)]));
    C_n = double(subs(C,[theta1,theta2,theta1_dot,theta2_dot],[y(i,1), y(i,2),y(i,3),y(i,4)]));
    G_n = double(subs(G,[theta1,theta2,theta1_dot,theta2_dot],[y(i,1), y(i,2),y(i,3),y(i,4)]));
    U = double(M_n*(-K *([y(i,1) ; y(i,2) ; y(i,3) ;y(i,4)]-[theta1_desired(i);theta2_desired(i);theta1_dot_desired(i);theta2_dot_desired(i)])+qd_ddot_sub)+C_n*[y(i,3) ;y(i,4)]+G_n);
%     disp("M")
%     disp(M)
%     disp("C")
%     disp(C)
%     disp("G")
%     disp(G)
    u1(i) = U(1);
    u2(i) = U(2);
    
end

%% Plotting 
fprintf("------------Plotting----------------\n");

figure
subplot(2,3,1)      
k1 = [y(:,1) ,theta1_desired];
plot(t,k1,'linewidth',2); 
title('Time vs Theta 1 and Theta 1 desired')
xlabel("Time in Seconds")
ylabel("Angle in Radians")
lgd = legend('theta1','theta1-desired');
lgd.FontSize = 10;

subplot(2,3,2)       
k2 = [y(:,2) ,theta2_desired];

plot(t,k2,'linewidth',2); 
title('Time vs Theta 2 and Theta 2 desired')
xlabel("Time in Seconds")
ylabel("Angle in Radians")
lgd = legend('theta2','theta2-desired');
lgd.FontSize = 10;

subplot(2,3,3)       
k3 = [y(:,3) ,theta1_dot_desired];

plot(t,k3,'linewidth',2); 
title('Time vs Theta 1-dot and Theta 1-dot desired')
xlabel("Time in Seconds")
ylabel("Angle in Radians")
lgd = legend('theta1-dot','theta1-dot-desired');
lgd.FontSize = 10;

subplot(2,3,4)       
k4 = [y(:,4) ,theta2_dot_desired];

plot(t,k4,'linewidth',2); 
title('Time vs Theta 2-dot and Theta 2-dot desired')
xlabel("Time in Seconds")
ylabel("Angle in Radians")
lgd = legend('theta2-dot','theta2-dot-desired');
lgd.FontSize = 10;

subplot(2,3,5)       
plot(t,u1,'linewidth',2); 
title('Time vs Torque at joint-1')
xlabel("Time in Seconds")
ylabel("Torque in N-m")
lgd = legend('Torque-1');
lgd.FontSize = 10;

subplot(2,3,6)       
plot(t,u2,'linewidth',2); 
title('Time vs Torque at joint-2')
xlabel("Time in Seconds")
ylabel("Torque in N-m")
lgd = legend('Torque-2');
lgd.FontSize = 10;


%% ODE Function 

function dz  = RR(t,z,K,M,C,G,m1,m2,l1,l2,r1,r2, I1 , I2 ,g,qd,qd_dot,qd_ddot)

dz = zeros(4,1);
z = num2cell(z);

[theta1 , theta2 , theta1_dot , theta2_dot] = deal(z{:}) ; 
qd= double(subs(qd));
qd_dot = double(subs(qd_dot));
qd_ddot = double(subs(qd_ddot));

if abs(theta1) > 2*pi
    theta1 = mod(theta1, 2*pi);
end 

if abs(theta2) > 2*pi
    theta2 = mod(theta2, 2*pi);
    
end 

% U = double(subs(M*(-K *[theta1 ; theta2 ; theta1_dot ;theta2_dot])+C*[theta1_dot;theta2_dot]+G)); %Control law with feedback linearization

U =  double(subs(M*(-K *([theta1 ; theta2 ; theta1_dot ;theta2_dot]-[qd;qd_dot])+qd_ddot)+C*[theta1_dot;theta2_dot]+G)); %Control law with trajectory tracking 
u1 = U(1,:);
u2 = U(2,:); 
dz(1) = theta1_dot; 
dz(2) = theta2_dot; 
dz(3) = (I2*u1 - I2*u2 + m2*r2^2*u1 - m2*r2^2*u2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) + g*l1*m2^2*r2^2*sin(theta1) + I2*g*l1*m2*sin(theta1) + I2*g*m1*r1*sin(theta1) - l1*m2*r2*u2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*m1*m2*r1*r2^2*sin(theta1) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);
dz(4) = -(I2*u1 - I1*u2 - I2*u2 - l1^2*m2*u2 - m1*r1^2*u2 + m2*r2^2*u1 - m2*r2^2*u2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1^3*m2^2*r2*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) - g*l1^2*m2^2*r2*sin(theta1 + theta2) - I1*g*m2*r2*sin(theta1 + theta2) + g*l1*m2^2*r2^2*sin(theta1) + I2*g*l1*m2*sin(theta1) + I2*g*m1*r1*sin(theta1) + l1*m2*r2*u1*cos(theta2) - 2*l1*m2*r2*u2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + 2*l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) + l1^2*m2^2*r2^2*theta2_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + g*l1^2*m2^2*r2*cos(theta2)*sin(theta1) - g*m1*m2*r1^2*r2*sin(theta1 + theta2) + I1*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*m1*m2*r1*r2^2*sin(theta1) + 2*l1^2*m2^2*r2^2*theta1_dot*theta2_dot*cos(theta2)*sin(theta2) + l1*m1*m2*r1^2*r2*theta1_dot^2*sin(theta2) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2) + g*l1*m1*m2*r1*r2*cos(theta2)*sin(theta1))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);
end

