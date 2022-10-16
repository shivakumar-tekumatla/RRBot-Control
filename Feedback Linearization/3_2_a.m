syms a0 a1 a2 a3 t 

q = a0+ a1*t + a2*t^2 + a3*t^3;
q_dot = jacobian(q,t);
t0=0;
tf =10;
eq1 = subs(q,[t],[t0]);
eq2 = subs(q,[t],[tf]);
eq3 = subs(q_dot,[t],[t0]);
eq4 = subs(q_dot,[t],[tf]);

sol = solve([eq1==180,eq2==0,eq3==0,eq4==0],[a0,a1,a2,a3]);
disp("a0")
sol.a0
disp("a1")
sol.a1
disp("a2")
sol.a2
disp("a3")
sol.a3