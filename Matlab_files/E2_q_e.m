function [E2] = E2_q_e(q)

  E2(1,1)=lf*sin(q(1) + q(2) - pi/2) - lf*sin(q(1) + q(3) - pi/2) + lt*sin(q(1) + q(2) + q(4) -...
          pi/2) - lt*cos(q(1) + q(3) + q(5));
  E2(1,2)=lf*sin(q(1) + q(2) - pi/2) + lt*sin(q(1) + q(2) + q(4) - pi/2);
  E2(1,3)=- lf*sin(q(1) + q(3) - pi/2) - lt*cos(q(1) + q(3) + q(5));
  E2(1,4)=lt*sin(q(1) + q(2) + q(4) - pi/2);
  E2(1,5)=-lt*cos(q(1) + q(3) + q(5));
  E2(1,6)=1;
  E2(1,7)=0;
  E2(2,1)=lf*cos(q(1) + q(2) - pi/2) - lf*cos(q(1) + q(3) - pi/2) + lt*cos(q(1) + q(2) + q(4) -...
          pi/2) + lt*sin(q(1) + q(3) + q(5));
  E2(2,2)=lf*cos(q(1) + q(2) - pi/2) + lt*cos(q(1) + q(2) + q(4) - pi/2);
  E2(2,3)=lt*sin(q(1) + q(3) + q(5)) - lf*cos(q(1) + q(3) - pi/2);
  E2(2,4)=lt*cos(q(1) + q(2) + q(4) - pi/2);
  E2(2,5)=lt*sin(q(1) + q(3) + q(5));
  E2(2,6)=0;
  E2(2,7)=1;

 