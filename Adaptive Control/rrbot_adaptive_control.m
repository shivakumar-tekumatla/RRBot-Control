clc;clear;close all;
%% Generalized coordinates 
% z =[theta1;theta2]
n = 2; %Number of generaized coordinates;
%% RR Bot control using Adaptive control law 
fprintf("------------Robust Law Implementation in Gazebo----------------\n");
i_linear_parametric_form;
j_adaptive_control_law;

rosshutdown;

rosinit;
j1_effort = rospublisher('/rrbot/joint1_effort_controller/command');%,'std_msgs/Float64');
j2_effort = rospublisher('/rrbot/joint2_effort_controller/command');%,'std_msgs/Float64');
JointStates = rossubscriber('/rrbot/joint_states');
tau1 = rosmessage(j1_effort);
tau2 = rosmessage(j2_effort);

tau1.Data = 0;
tau2.Data = 0;
send(j1_effort,tau1);
send(j2_effort,tau2);
client = rossvcclient('/gazebo/set_model_configuration');
req = rosmessage(client);
req.ModelName = 'rrbot';
req.UrdfParamName = 'robot_description';
req.JointNames = {'joint1','joint2'};


t=0;
time = [];
joint1_pos = [];
joint2_pos = [];
joint1_vel = [];
joint2_vel = [];
input1     = [];
input2     = [];
estimated_values=[];
desired_trajectory=[];
actual_trajectory=[];
errors =[];
req.JointPositions = [deg2rad(180), deg2rad(90)];
resp = call(client,req,'Timeout',3);
v=[0;0];
rho_0=0.001;
phi=0.01;
tic;
while (t<10)
    t = toc;
    time = [time,t];
    theta1_desired =(pi*t^3)/500 - (3*pi*t^2)/100 + pi;
    theta2_desired =  (pi*t^3)/1000 - (3*pi*t^2)/200 + pi/2;
    theta1_dot_desired = (3*pi*t^2)/500 - (3*pi*t)/50;
    theta2_dot_desired =    (3*pi*t^2)/1000 - (3*pi*t)/100;
    xref = [theta1_desired;theta2_desired;theta1_dot_desired;theta2_dot_desired];
    desired_trajectory = [desired_trajectory,xref];
    % read the joint states
    jointData = receive(JointStates);
    x= [jointData.Position;jointData.Velocity];
    actual_trajectory=[actual_trajectory,x];
    
    q1 = x(1);
    q2 = x(2);
    dq1 = x(3);
    dq2 = x(4);
    e = x-xref;
    errors=[errors,e];

    ddqd = [(3*pi*t)/250 - (3*pi)/50;
                (3*pi*t)/500 - (3*pi)/100];
    
    K = [Kp Kd];
    
    v = ddqd- K*e ;
    ddq1 = v(1);
    ddq2 = v(2);
    a_hat  = alpha_hat(1); 
    b_hat  = alpha_hat(2); 
    d_hat  = alpha_hat(3); 

    
    Mmat_hat= [a_hat+2*b_hat*cos(q2), d_hat+b_hat*cos(q2); d_hat+b_hat*cos(q2), d_hat];
 
    Y =[ddq1, - sin(q2)*dq2^2 - 2*dq1*sin(q2)*dq2 + cos(q2)*(2*ddq1 + ddq2),        ddq2, -(981*sin(q1))/100, -(981*sin(q1 + q2))/100;
        0,                                  sin(q2)*dq1^2 + ddq1*cos(q2), ddq1 + ddq2,                  0, -(981*sin(q1 + q2))/100];
    
    [t_e,alpha_hat] = ode45(@(t,x_alpha) esimates(t,x_alpha,Y,Mmat_hat,B,P,e,gamma), [time(end)-0.0001,t],[alpha_hat]);
    
    alpha_hat=alpha_hat(end,:)';
    estimated_values = [estimated_values,alpha_hat];
    tau= Y*alpha_hat; 

    u1 =tau(1);
    u2 =tau(2);
    input1=[input1,u1];
    input2=[input2,u2];
    tau1.Data = u1;%input1;
    tau2.Data = u2;%5;%input2;S
    send(j1_effort,tau1);
    send(j2_effort,tau2);
     
    
end
tau1.Data = 0;
tau2.Data = 0;
send(j1_effort,tau1);
send(j2_effort,tau2);% disconnectfromroscore

rosshutdown;

%%
% Plotting Joint Variables
fprintf("------------Plotting Trajectories and Control Inputs----------------\n");
disp("Plotting Results")


figure(1);
plot(time,actual_trajectory,'linewidth',2);
title('Time vs Theta1, Theta2, Theta1-dot, Theta2-dot')
lgd = legend('Theta1', 'Theta2', 'Theta1-dot', 'Theta2-dot');
lgd.FontSize = 14;

figure(2);
subplot(2,2,1)      
plot(time',actual_trajectory(1,:)',"b-",time',desired_trajectory(1,:)',"r--",'linewidth',2);   
lgd = legend('Theta-1','Theta1-desired');
lgd.FontSize = 14;
title('Time vs Theta-1')
xlabel("Time in Seconds")
ylabel("Angle in radians")

subplot(2,2,2)      
plot(time',actual_trajectory(2,:)',"b-",time',desired_trajectory(2,:)',"r--",'linewidth',2);   
lgd = legend('Theta-2','Theta2-desired');
lgd.FontSize = 14;
title('Time vs Theta-2')
xlabel("Time in Seconds")
ylabel("Angle in radians")

subplot(2,2,3)      
plot(time',actual_trajectory(3,:)',"b-",time',desired_trajectory(3,:)',"r--",'linewidth',2);      
lgd = legend('Theta-1-dot','Theta1-dot-desired');
lgd.FontSize = 14;
title('Time vs Theta-1-dot')
xlabel("Time in Seconds")
ylabel("Velocity in radians per second")

subplot(2,2,4)      
plot(time',actual_trajectory(4,:)',"b-",time',desired_trajectory(4,:)',"r--",'linewidth',2);   
lgd = legend('Theta-2-dot','Theta2-dot-desired');
lgd.FontSize = 14;
title('Time vs Theta-2-dot')
xlabel("Time in Seconds")
ylabel("Velocity in radians per second")

figure(3);
subplot(2,1,1);
plot(time,input1,'linewidth',2);
title("Joint-1 Torque")
lgd = legend("Joint-1 Torque");
lgd.FontSize = 14;
xlabel("Time in Seconds")
ylabel("Torque in N-m")
subplot(2,1,2);
plot(time,input2,'linewidth',2);
title("Joint-2 Torque")
lgd = legend("Joint-2 Torque");
lgd.FontSize = 14;
xlabel("Time in Seconds")
ylabel("Torque in N-m")

figure(4)
plot(time,errors,'linewidth',2);
title('Time vs Errors')
lgd = legend('Theta1-Error', 'Theta2-Error', 'Theta1-dot-Error', 'Theta2-dot-Error');
lgd.FontSize = 14;

figure(5)
plot(time,estimated_values,'LineWidth',2);
title("Time vs Estimates")
lgd = legend("alpha-1","alpha-2","alpha-3","alpha-4","alpha-5");
ldg.Fontsize =14;


%% ODE 45 to find updated estimates 
function dx = esimates(t,alpha_hat,Y,Mmat_hat,B,P,e,gamma)
        
    PHI = inv(Mmat_hat)*Y;
%     disp(PHI);
    
    alpha_hat_dot = -inv(gamma)*(PHI'*B'*P*e);
    
    dx(1,1) = alpha_hat_dot(1);
    dx(2,1) = alpha_hat_dot(2);
    dx(3,1) = alpha_hat_dot(3);
    dx(4,1) = alpha_hat_dot(4);
    dx(5,1) = alpha_hat_dot(5);

end