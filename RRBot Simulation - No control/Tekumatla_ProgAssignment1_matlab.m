%{
generalized coordinates theta1 , theta2 ===>  q = transpose[theta1 theta2]

generailzed forces are u1, u2 ====> u = transpose[u1 u2]

There are two different links, hence we calculate energies for both of
them separately 

KE1 is kinetic energy of link 1 

KE2 is kinetic energy of link 2 

PE1 is potential energy if link1 

PE2 is potential energy of link2 

L = (KE1+KE2) - (PE1+PE2) 

%}

%Problem 1.2 a 

syms l1 l2 r1 r2 theta1 theta2 I1 I2 m1 m2 theta1_dot theta2_dot g theta1_ddot theta2_ddot u1 u2 

KE1 = 0.5*m1*(r1*theta1_dot)^2 + 0.5*I1*theta1_dot^2;

PE1 = m1*g*r1*cos(theta1);

KE2 = 0.5*m2*(l1*theta1_dot)^2 + 0.5*m2*(r2*(theta1_dot+theta2_dot))^2 + m2*l1*theta1_dot*r2*(theta1_dot+theta2_dot)*cos(theta2) + 0.5*I2*(theta1_dot+theta2_dot)^2;

PE2 = m2*g*l1*cos(theta1) + m2*g*r2*cos(theta1+theta2);

L = (KE1 + KE2) - (PE1+PE2);
disp("Legrangian is")
disp(L);

partial_diff =  jacobian(L , [theta1,theta1_dot,theta2,theta2_dot]);

partial_diff(2) = jacobian(partial_diff(2),[theta1;theta2;theta1_dot;theta2_dot])*[theta1_dot;theta2_dot;theta1_ddot;theta2_ddot];

partial_diff(4) = jacobian(partial_diff(4),[theta1;theta2;theta1_dot;theta2_dot])*[theta1_dot;theta2_dot;theta1_ddot;theta2_ddot];

eq1 = partial_diff(2)-partial_diff(1)-u1;

eq2 = partial_diff(4)-partial_diff(3)-u2;

disp("Dynamic Equations are")

disp("Equation1 ")
disp(eq1);

disp("Equation2 ")
disp(eq2);


% Problem 1.2 b -- State Space Model 

X = sym('X' , [4,1]);
X(1) = theta1;
X(2) = theta2;
X(3) =  theta1_dot;
X(4) = theta2_dot;
sol = solve( [eq1==0 , eq2==0] , [theta1_ddot,theta2_ddot]);

disp("Solutions are");

disp("theta1_ddot") ;
disp(sol.theta1_ddot);
disp("theta2_ddot");
disp(sol.theta2_ddot);


% writing in manipulator form 
% M*q_ddot+C*q_dot+g = F


M =  [m2*l1^2+2*m2*l2*r2*cos(theta2)+m1*r1^2+m2*r2^2+I1+I2  m2*r2^2+m2*r2*l1*cos(theta2)+I2;
    m2*r2^2+m2*l1*r2*cos(theta2)+I2                         m2*r2^2+I2];

C = [0 l1*m2*r2*theta1_dot*sin(theta2)+l1*m2*r2*sin(theta2)*(theta1_dot+theta2_dot);
    l1*m2*r2*(theta1_dot+theta2_dot)*sin(theta2)-l1*m2*r2*theta2_dot*sin(theta2) 0];

g = [-g*l1*m2*sin(theta1)-g*m1*r1*sin(theta1)-m2*g*r2*sin(theta1+theta2);
    -g*m2*r2*sin(theta1+theta2)];

F = [u1;u2];

MT = transpose(M);

if MT==M 
    
    disp("M is a symmetric matrix")
else 
    disp("M is NOT  a symmetric matrix")
end

m1=1 ; m2= 1;  l1=1 ; l2= 1; r1=0.45 ; r2= 0.45;
I1=0.084; I2= 0.084;  g= 9.81; % Used same manipulator parameters as given in the problem

theta2 = pi/4 ; % Enter any theta_2 value here

M = subs(M);

eig_val = eig(M);

if isAlways(eig_val>0)
    
     disp("M is a positive definite matrix")

else
    disp("M is NOT a positive definite matrix")
    
end 
    
% Problem 1.2 C --- Simulation  

figure(1);
[t,y] = ode45(@RR , [0,10],[deg2rad(200),deg2rad(125),0,0]);

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


function dz  = RR(t,z)

m1=1 ; m2= 1;  l1=1 ; l2= 1; r1=0.45 ; r2= 0.45;
I1=0.084; I2= 0.084;  g= 9.81;

dz = zeros(4,1);
z = num2cell(z);

[theta1 , theta2 , theta1_dot , theta2_dot] = deal(z{:}) ; 

if abs(theta1) > 2*pi
    theta1 = mod(theta1, 2*pi);
end 

if abs(theta2) > 2*pi
    theta2 = mod(theta2, 2*pi);
    
end 

u1 = 0; u2 =0; 

dz(1) = theta1_dot; 

dz(2) = theta2_dot; 

dz(3) = (I2*u1 - I2*u2 + m2*r2^2*u1 - m2*r2^2*u2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) + g*l1*m2^2*r2^2*sin(theta1) + I2*g*l1*m2*sin(theta1) + I2*g*m1*r1*sin(theta1) - l1*m2*r2*u2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*m1*m2*r1*r2^2*sin(theta1) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);

dz(4) = -(I2*u1 - I1*u2 - I2*u2 - l1^2*m2*u2 - m1*r1^2*u2 + m2*r2^2*u1 - m2*r2^2*u2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1^3*m2^2*r2*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) - g*l1^2*m2^2*r2*sin(theta1 + theta2) - I1*g*m2*r2*sin(theta1 + theta2) + g*l1*m2^2*r2^2*sin(theta1) + I2*g*l1*m2*sin(theta1) + I2*g*m1*r1*sin(theta1) + l1*m2*r2*u1*cos(theta2) - 2*l1*m2*r2*u2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + 2*l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) + l1^2*m2^2*r2^2*theta2_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + g*l1^2*m2^2*r2*cos(theta2)*sin(theta1) - g*m1*m2*r1^2*r2*sin(theta1 + theta2) + I1*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*m1*m2*r1*r2^2*sin(theta1) + 2*l1^2*m2^2*r2^2*theta1_dot*theta2_dot*cos(theta2)*sin(theta2) + l1*m1*m2*r1^2*r2*theta1_dot^2*sin(theta2) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2) + g*l1*m1*m2*r1*r2*cos(theta2)*sin(theta1))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);
end