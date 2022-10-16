 % g) RR Bot control using robust control law 
clc;clear;close all;
clear classes;
rehash toolboxcache;
fprintf("------------Robust Law Implementation in Gazebo----------------\n");
c_robust_control_law;
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
desired_trajectory=[];
actual_trajectory=[];
errors =[];
req.JointPositions = [deg2rad(200), deg2rad(125)];
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
%      e1(i) = e(1);
%      e2(i) = e(2);
%      e1_dot(i) = e(3);
%      e2_dot(i) = e(4);
     
    errors=[errors,e];
     
    Mmat  = [(9*cos(q2))/10 + 1573/1000  (9*cos(q2))/20 + 573/2000;
             (9*cos(q2))/20 + 573/2000                   573/2000];
    Cmat =[-(9*dq2*sin(q2))/20, -(9*sin(q2)*(dq1 + dq2))/20;
             (9*dq1*sin(q2))/20,                           0];
    Gmat =[- (8829*sin(q1 + q2))/2000 - (28449*sin(q1))/2000;
            -(8829*sin(q1 + q2))/2000];
    %Nominal parameters
    
    Mmat_hat =[(27*cos(q2))/40 + 4719/4000, (27*cos(q2))/80 + 1719/8000;
        (27*cos(q2))/80 + 1719/8000,                   1719/8000];
    Cmat_hat =[-(27*dq2*sin(q2))/80, -(27*sin(q2)*(dq1 + dq2))/80;
            (27*dq1*sin(q2))/80,                            0];
    Gmat_hat =[- (26487*sin(q1 + q2))/8000 - (85347*sin(q1))/8000;
                        -(26487*sin(q1 + q2))/8000];
     ddqd = [(3*pi*t)/250 - (3*pi)/50;
            (3*pi*t)/500 - (3*pi)/100];
    Cmat_tilde = Cmat_hat-Cmat;
    Gmat_tilde = Gmat_hat-Gmat;
    
    uncertainity = (inv(Mmat)*Mmat_hat-eye(size(Mmat)))*v+inv(Mmat)*(Cmat_tilde*[dq1;dq2]+Gmat_tilde);
%     disp(uncertainity);
    rho = norm(uncertainity)+rho_0; 
     if phi>0
        if norm(e'*P*B)>phi
            vr = -rho*(e'*P*B)/(norm(e'*P*B));
        else
            vr = -rho*(e'*P*B)/phi;
        end
    else
        if norm(e'*P*B)~=0
            vr = -rho*(e'*P*B)/(norm(e'*P*B));
        else
            vr=0;
        end
    end
     
     v = ddqd- K*e +vr';
%     disp(v)
    tau = Mmat_hat*v + Cmat_hat*[dq1;dq2]+Gmat_hat;
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