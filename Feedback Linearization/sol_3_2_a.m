clc;clear;close all;
syms a0 a1 a2 a3 b0 b1 b2 b3 t
%% Solving Cubic polynomials
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

sol = solve([eq1==180,eq2==0,eq3==0,eq4==0,eq5==90,eq6==0,eq7==0,eq8==0],[a0,a1,a2,a3,b0,b1,b2,b3]);

sprintf("a0 = %f",sol.a0)
sprintf("a1 = %f",sol.a1)
sprintf("a2 = %f",sol.a2)
sprintf("a3 = %f",sol.a3)
sprintf("b0 = %f",sol.b0)
sprintf("b1 = %f",sol.b1)
sprintf("b2 = %f",sol.b2)
sprintf("b3 = %f",sol.b3)

cubic_1 = subs(q1,[a0,a1,a2,a3],[sol.a0,sol.a1,sol.a2,sol.a3]);
cubic_2 = subs(q2,[b0,b1,b2,b3],[sol.b0,sol.b1,sol.b2,sol.b3]);
disp("Cubic Polynomial for 1st Joint")
disp(cubic_1)
disp("Cubic Polynomial for 2nd Joint")
disp(cubic_2)