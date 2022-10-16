
%state space model of the 2R manipulator system 

syms theta1 theta2  theta1_dot theta2_dot  theta1_ddot theta2_ddot u1 u2 
syms x1 x2 x3 x4 k11 k12 k13 k14 k21 k22 k23 k24

%Equation Motions taken from the Programming assignment 1 

%2.2(a) 
m1 = 1; m2= 1; l1=1;l2=1;r1=0.45;r2=0.45; I1 = 0.084; I2 = 0.084; g =9.81;

eq1 = theta1_ddot*(m2*l1^2 + 2*m2*cos(theta2)*l1*r2 + m1*r1^2 + m2*r2^2 + I1 + I2) - theta2_dot*(l1*m2*r2*theta1_dot*sin(theta2) + l1*m2*r2*sin(theta2)*(theta1_dot + theta2_dot)) - u1 + theta2_ddot*(m2*r2^2 + l1*m2*cos(theta2)*r2 + I2) - g*m2*r2*sin(theta1 + theta2) - g*l1*m2*sin(theta1) - g*m1*r1*sin(theta1);

eq2 = theta2_ddot*(m2*r2^2 + I2) - u2 + theta1_ddot*(m2*r2^2 + l1*m2*cos(theta2)*r2 + I2) - g*m2*r2*sin(theta1 + theta2) + l1*m2*r2*theta1_dot*sin(theta2)*(theta1_dot + theta2_dot) - l1*m2*r2*theta1_dot*theta2_dot*sin(theta2);

eq_eq1 = subs(eq1,[theta1_dot,theta2_dot,theta1_ddot,theta2_ddot,u1,u2],[0,0,0,0,0,0]) ; %subs(eq1, [theta1_dot,theta2_dot,theta1_ddot,theta2_ddot,u1,u2,m1,m2,l1,l2,r1,r2,I1,I2,g],[0,0,0,0,0,0,1,1,1,1,0.45,0.45,0.084,0.084,9.81]);
eq_eq2 = subs(eq2,[theta1_dot,theta2_dot,theta1_ddot,theta2_ddot,u1,u2],[0,0,0,0,0,0]);  %subs(eq2, [theta1_dot,theta2_dot,theta1_ddot,theta2_ddot,u1,u2,m1,m2,l1,l2,r1,r2,I1,I2,g],[0,0,0,0,0,0,1,1,1,1,0.45,0.45,0.084,0.084,9.81]);

sol_eq = solve(eq_eq1==0, eq_eq2==0,[theta1,theta2]);

eq_points = [sol_eq.theta1,sol_eq.theta2,[0;0;0],[0;0;0],[0;0;0],[0;0;0]]; %[theta1,theta2,theta1_dot,theta2_dot,u1,u2]
disp("Equilibriam points are");
disp(eq_points);
%2.2 (b)  Jacobian Linearization 

% statespace model 

X_dot = sym('X_dot' , [4,1]);
X_dot(1) = theta1_dot;
X_dot(2) = theta2_dot;
X_dot(3) =  (I2*u1 - I2*u2 + m2*r2^2*u1 - m2*r2^2*u2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) + g*l1*m2^2*r2^2*sin(theta1) + I2*g*l1*m2*sin(theta1) + I2*g*m1*r1*sin(theta1) - l1*m2*r2*u2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*m1*m2*r1*r2^2*sin(theta1) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);
X_dot(4) = -(I2*u1 - I1*u2 - I2*u2 - l1^2*m2*u2 - m1*r1^2*u2 + m2*r2^2*u1 - m2*r2^2*u2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1^3*m2^2*r2*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) - g*l1^2*m2^2*r2*sin(theta1 + theta2) - I1*g*m2*r2*sin(theta1 + theta2) + g*l1*m2^2*r2^2*sin(theta1) + I2*g*l1*m2*sin(theta1) + I2*g*m1*r1*sin(theta1) + l1*m2*r2*u1*cos(theta2) - 2*l1*m2*r2*u2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + 2*l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) + l1^2*m2^2*r2^2*theta2_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + g*l1^2*m2^2*r2*cos(theta2)*sin(theta1) - g*m1*m2*r1^2*r2*sin(theta1 + theta2) + I1*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*m1*m2*r1*r2^2*sin(theta1) + 2*l1^2*m2^2*r2^2*theta1_dot*theta2_dot*cos(theta2)*sin(theta2) + l1*m1*m2*r1^2*r2*theta1_dot^2*sin(theta2) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2) + g*l1*m1*m2*r1*r2*cos(theta2)*sin(theta1))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);

A = jacobian(X_dot,[theta1,theta2,theta1_dot,theta2_dot]);

B = jacobian(X_dot, [u1,u2]);

% Linearization around the equilibrium points 

for eq_point = eq_points'
    disp('Equilibrium Point(Theta_1,Theta_2,Theta_1_dot,Theta_2_dot,u1,u2');
    disp(eq_point');
    A_eq = double(subs(A,[theta1,theta2,theta1_dot,theta2_dot,u1,u2],eq_point')); 
    B_eq = double(subs(B,[theta1,theta2,theta1_dot,theta2_dot,u1,u2],eq_point')); 
    
    disp('State Matrix');
    disp(A_eq);
%     disp('Eigen Values');
    
    %Stability check 
    eig_vals = eig(A_eq);
%     disp(eig_vals);
    count =0;
    for i= eig_vals'
%         disp(real(i))
        if real(i)>0
            count= count+1;
%             disp(count)
        end 
    end 
    disp('Input Matrix');
    disp(B_eq);
%     stable = all(real(eig_vals)<=0);
%     disp(real(eig_vals))

    if count == 0
        disp('The system is stable at the above equilibrium point')
    else
        disp('The system is NOT stable at the above equilibrium point')
    end 
    
    %check controllability 
    if rank(ctrb(A_eq,B_eq))==4
        disp("The system is  Controllable at the above equilibrium point")
    else
        disp("The system is  NOT Controllable at the above equilibrium point")
    end 
    disp("_____________________________________________________________________________________")
end 

%Upward Equilibrium point is [ 0,0,0,0] 
%state Feedback at this equilibrium and eigen values at [-2-i,-2+i,-3,-5]

A_up = double(subs(A,[theta1,theta2,theta1_dot,theta2_dot],[0,0,0,0]));

B_up = double(subs(B,[theta1,theta2,theta1_dot,theta2_dot],[0,0,0,0]));

lambda = [-2-i,-2+i,-1,-0.5]; %Eigen Value placement 

Kn = place(A_up,B_up,lambda);

%observer 

C= [1,1,0,0];

if rank(obsv(A_up,C))==4
    
    disp("System is Observable");
else
    disp("No Observable")
    
end 

L = place(A_up',C',lambda);
figure(1);

% [t,y] = ode45(@RR , [0,10],[pi/6,pi/4,0,0],Kn);

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

function dz  = RR(t,z,Kn,m1,m2,l1,l2,r1,r2,I1,I2,g)

% m1=1 ; m2= 1;  l1=1 ; l2= 1; r1=0.45 ; r2= 0.45;
% I1=0.084; I2= 0.084;  g= 9.81;

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