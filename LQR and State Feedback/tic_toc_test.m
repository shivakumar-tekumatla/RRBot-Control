
i= 1;
while i <= size(input1)
    
    % read the joint states
    jointData = receive(JointStates);
    % inspect the "jointData" variableinMATLAB to get familiar with itsstructure
    % design your state feedback controller in the following
    tau1.Data = input1(i);
    tau2.Data = input2(i);
    send(j1_effort,tau1);
    send(j2_effort,tau2);
    % you can sample data here to plot at the end
    i = i+1 ;
end