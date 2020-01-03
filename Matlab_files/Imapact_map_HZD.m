function y= Impact_map_HZD(x)
qsd=[x(6);x(7);x(8);x(9);x(10)];
De=D_q_e(x);
E2=E2_q_e(x);
A=[De -E2';E2 zeros(2,2)];
B = [De*[qsd;zeros(2,1)]; zeros(2,1)];
q_plus = A\B;
y = [x(1); x(3); x(2); x(5); x(4); q_plus(1); q_plus(3); q_plus(2); q_plus(5); q_plus(4)];