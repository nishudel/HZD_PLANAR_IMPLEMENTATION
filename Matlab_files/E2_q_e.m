function [E2] = E2_q_e(q)

  E2(1,1)=(2*sin(q(1) + q(2) - pi/2))/5 - (2*cos(q(1) + q(3) + q(5)))/5 - (2*sin(q(1) + q(3) -...
          pi/2))/5 + (2*sin(q(1) + q(2) + q(4) - pi/2))/5;
  E2(1,2)=(2*sin(q(1) + q(2) - pi/2))/5 + (2*sin(q(1) + q(2) + q(4) - pi/2))/5;
  E2(1,3)=- (2*cos(q(1) + q(3) + q(5)))/5 - (2*sin(q(1) + q(3) - pi/2))/5;
  E2(1,4)=(2*sin(q(1) + q(2) + q(4) - pi/2))/5;
  E2(1,5)=-(2*cos(q(1) + q(3) + q(5)))/5;
  E2(1,6)=1;
  E2(1,7)=0;
  E2(2,1)=(2*sin(q(1) + q(3) + q(5)))/5 + (2*cos(q(1) + q(2) - pi/2))/5 - (2*cos(q(1) + q(3) -...
          pi/2))/5 + (2*cos(q(1) + q(2) + q(4) - pi/2))/5;
  E2(2,2)=(2*cos(q(1) + q(2) - pi/2))/5 + (2*cos(q(1) + q(2) + q(4) - pi/2))/5;
  E2(2,3)=(2*sin(q(1) + q(3) + q(5)))/5 - (2*cos(q(1) + q(3) - pi/2))/5;
  E2(2,4)=(2*cos(q(1) + q(2) + q(4) - pi/2))/5;
  E2(2,5)=(2*sin(q(1) + q(3) + q(5)))/5;
  E2(2,6)=0;
  E2(2,7)=1;

 