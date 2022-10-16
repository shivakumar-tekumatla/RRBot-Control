 % g) RR Bot control using Gazebo 

fprintf("------------Feedback Control Law Implementation in Gazebo----------------\n");

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
req.JointPositions = [deg2rad(30), deg2rad(45)];
resp = call(client,req,'Timeout',3);

tic;
t=0;
time = [];
joint1_pos = [];
joint2_pos = [];
joint1_vel = [];
joint2_vel = [];
input1     = [];
input2     = [];

while (t<10)
    t = toc;
    time = [time,t];

    % read the joint states
    jointData = receive(JointStates);
    % inspect the "jointData" variable in MATLAB to get familiar with itsstructure
    % design your state feedback controller in the following
    u1 = -Kn(1,:) *[jointData.Position;jointData.Velocity];
    u2 = -Kn(2,:) *[jointData.Position;jointData.Velocity];
    tau1.Data = u1;%input1;
    tau2.Data = u2;%5;%input2;S
    send(j1_effort,tau1);
    send(j2_effort,tau2);
    showdetails(jointData);
    % you can sample data here to plot at the end
    joint1_pos=[joint1_pos,jointData.Position(1)];
    joint2_pos=[joint2_pos,jointData.Position(2)];
    joint1_vel=[joint1_vel,jointData.Velocity(1)];
    joint2_vel=[joint2_vel,jointData.Velocity(2)];
    input1    =[input1,u1];
    input2    =[input2,u2];
    
%     pause(10/177,t);
    
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

figure
subplot(2,2,1)      
plot(time,joint1_pos,'linewidth',2);          
title('Time vs Theta 1')
xlabel("Time in Seconds")
ylabel("Angle in Radians")

subplot(2,2,2)       
plot(time,joint2_pos,'linewidth',2);       
title('Time vs Theta 2')
xlabel("Time in Seconds")
ylabel("Angle in Radians")

subplot(2,2,3)       
plot(time,joint1_vel,'linewidth',2);      
title('Time vs Theta 1-Dot')
xlabel("Time in Seconds")
ylabel("Angular Velocity in Radians per second ")

subplot(2,2,4)       
plot(time,joint2_vel,'linewidth',2);        
title('Time vs Theta 2-Dot')
xlabel("Time in Seconds")
ylabel("Angular Velocity in Radians per second ")

%Plot Forces/Torques
figure
subplot(2,1,1)       
plot(time,input1,'linewidth',2);        
title('Time vs Joint -1 Torque')
xlabel("Time in Seconds")
ylabel("Torque(N-m) ")

subplot(2,1,2)       
plot(time,input2,'linewidth',2);        
title('Time vs Joint -2 Torque')
xlabel("Time in Seconds")
ylabel("Torque(N-m) ")
