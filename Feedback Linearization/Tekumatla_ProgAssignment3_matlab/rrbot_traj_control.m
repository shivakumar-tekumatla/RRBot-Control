clc;clear;close all;
% g) RR Bot control using Gazebo 
fprintf("------------Feedback Control Law Implementation in Gazebo----------------\n");

rosshutdown;

rosinit;
j1_effort = rospublisher('/rrbot/joint1_effort_controller/command');%,'std_msgs/Float64');
j2_effort = rospublisher('/rrbot/joint2_effort_controller/command');%,'std_msgs/Float64');
JointStates = rossubscriber('/rrbot/joint_states');
tau1 = rosmessage(j1_effort);
tau2 = rosmessage(j2_effort);
joints = rosmessage(JointStates);

tau1.Data = 0;
tau2.Data = 0;

send(j1_effort,tau1);
send(j2_effort,tau2);
client = rossvcclient('/gazebo/set_model_configuration');
req = rosmessage(client);
req.ModelName = 'rrbot';
req.UrdfParamName = 'robot_description';
req.JointNames = {'joint1','joint2'};
theta1 = deg2rad(200);
theta2 = deg2rad(125);
theta1_dot=0;
theta2_dot=0;
req.JointPositions = [theta1,theta2];
resp = call(client,req,'Timeout',3);

time = [];
joint1_pos = [];
joint2_pos = [];
joint1_vel = [];
joint2_vel = [];
input1     = [];
input2     = [];
% qd = [cubic_1;cubic_2];
% qd_dot = jacobian(qd);
% qd_ddot = jacobian(qd_dot);
m1 = 1; m2= 1; l1=1;l2=1;r1=0.45;r2=0.45; I1 = 0.084; I2 = 0.084; g =9.81;

M = [m2*l1^2+2*m2*cos(theta2)*l2*r2+m1*r1^2+m2*r2^2+I1+I2   m2*r2^2+m2*r2*l1*cos(theta2)+I2;
    m2*r2^2+m2*r2*l1*cos(theta2)+I2                          m2*r2^2+I2];

C = [0                                         l1*m2*r2*theta1_dot*sin(theta2)+l1*m2*r2*sin(theta2)*(theta1_dot+theta2_dot);
    l1*m2*r2*(theta1_dot+theta2_dot)*sin(theta2)-l1*m2*r2*theta2_dot*sin(theta2)            0];
G = [-g*l1*m2*sin(theta1)-g*m1*r1*sin(theta1)-m2*g*r2*sin(theta1+theta2);
    -g*m2*r2*sin(theta1+theta2)];
K =[8.4365    4.7619    5.7267    1.9061;
   -5.9615    8.1920   -1.9986    5.7733];

% K = [2 0.75 2 0.75;
%     4 1.5 4 1.5];


tic;
t=0;
while (t<10)
    t = toc;
    time = [time,t];
    
    jointData = receive(JointStates);
    theta1 = jointData.Position(1);
    theta2 = jointData.Position(2);
    theta1_dot = jointData.Velocity(1);
    theta2_dot = jointData.Velocity(2);
    
    if abs(theta1) > 2*pi
        theta1 = mod(theta1, 2*pi);
    end 

    if abs(theta2) > 2*pi
        theta2 = mod(theta2, 2*pi);
    end 
    
%     ddddd=[theta1;theta2]
    M = [m2*l1^2+2*m2*cos(theta2)*l2*r2+m1*r1^2+m2*r2^2+I1+I2   m2*r2^2+m2*r2*l1*cos(theta2)+I2;
    m2*r2^2+m2*r2*l1*cos(theta2)+I2                          m2*r2^2+I2];

    C = [0                                         l1*m2*r2*theta1_dot*sin(theta2)+l1*m2*r2*sin(theta2)*(theta1_dot+theta2_dot);
        l1*m2*r2*(theta1_dot+theta2_dot)*sin(theta2)-l1*m2*r2*theta2_dot*sin(theta2)            0];
    G = [-g*l1*m2*sin(theta1)-g*m1*r1*sin(theta1)-m2*g*r2*sin(theta1+theta2);
        -g*m2*r2*sin(theta1+theta2)];
    
    z = [theta1;theta2;theta1_dot;theta2_dot];
    zd =[(pi*t^3)/500 - (3*pi*t^2)/100 + pi;(pi*t^3)/1000 - (3*pi*t^2)/200 + pi/2;(3*pi*t^2)/500 - (3*pi*t)/50;(3*pi*t^2)/1000 - (3*pi*t)/100];
    vd = [(3*pi*t)/250 - (3*pi)/50;(3*pi*t)/500 - (3*pi)/100];
    virtual_input = -K*(z-zd)+vd;
    U = double(subs(M*virtual_input+C*z(3:4)+G));
    u1 = U(1);
    u2 = U(2);
    tau1.Data = u1;%input1;
    tau2.Data = u2;%5;%input2;S
    send(j1_effort,tau1);
    send(j2_effort,tau2);
%     showdetails(jointData);
    joint1_pos=[joint1_pos,theta1];
    joint2_pos=[joint2_pos,theta2];
    joint1_vel=[joint1_vel,theta1_dot];
    joint2_vel=[joint2_vel,theta2_dot];
    input1    =[input1,u1];
    input2    =[input2,u2];
    
end
tau1.Data = 0;
tau2.Data = 0;
send(j1_effort,tau1);
send(j2_effort,tau2);% disconnectfromroscore

rosshutdown;

%% Desired Trajectories comuptation


theta1_desired =zeros(length(time),1);
theta2_desired =zeros(length(time),1);
theta1_dot_desired =zeros(length(time),1);
theta2_dot_desired =zeros(length(time),1);
syms t;

qd = [(pi*t^3)/500 - (3*pi*t^2)/100 + pi;
      (pi*t^3)/1000 - (3*pi*t^2)/200 + pi/2];
qd_dot = [(3*pi*t^2)/500 - (3*pi*t)/50;
          (3*pi*t^2)/1000 - (3*pi*t)/100];

for i=1:length(time)
    t= time(i);
    qd_sub = subs(qd,t);
    qd_dot_sub = subs(qd_dot,t);
    theta1_desired(i) = qd_sub(1);
    theta2_desired(i) = qd_sub(2);
    theta1_dot_desired(i) = qd_dot_sub(1);
    theta2_dot_desired(i) = qd_dot_sub(2); 
end


%%
% Plotting Joint Variables
fprintf("------------Plotting Trajectories and Control Inputs----------------\n");
disp("Plotting Results")


figure
subplot(2,3,1)      
k1 = [joint1_pos',theta1_desired];
plot(time,k1,'linewidth',2); 
title('Time vs Theta 1 and Theta 1 desired')
xlabel("Time in Seconds")
ylabel("Angle in Radians")
lgd = legend('theta1','theta1-desired');
lgd.FontSize = 10;

subplot(2,3,2)       
k2 = [joint2_pos',theta2_desired];

plot(time,k2,'linewidth',2); 
title('Time vs Theta 2 and Theta 2 desired')
xlabel("Time in Seconds")
ylabel("Angle in Radians")
lgd = legend('theta2','theta2-desired');
lgd.FontSize = 10;

subplot(2,3,3)       
k3 = [joint1_vel',theta1_dot_desired];

plot(time,k3,'linewidth',2); 
title('Time vs Theta 1-dot and Theta 1-dot desired')
xlabel("Time in Seconds")
ylabel("Angle in Radians")
lgd = legend('theta1-dot','theta1-dot-desired');
lgd.FontSize = 10;

subplot(2,3,4)       
k4 = [joint2_vel',theta2_dot_desired];

plot(time,k4,'linewidth',2); 
title('Time vs Theta 2-dot and Theta 2-dot desired')
xlabel("Time in Seconds")
ylabel("Angle in Radians")
lgd = legend('theta2-dot','theta2-dot-desired');
lgd.FontSize = 10;

subplot(2,3,5)       
plot(time,input1,'linewidth',2); 
title('Time vs Torque at joint-1')
xlabel("Time in Seconds")
ylabel("Torque in N-m")
lgd = legend('Torque-1');
lgd.FontSize = 10;

subplot(2,3,6)       
plot(time,input2,'linewidth',2); 
title('Time vs Torque at joint-2')
xlabel("Time in Seconds")
ylabel("Torque in N-m")
lgd = legend('Torque-2');
lgd.FontSize = 10;