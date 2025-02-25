function [D] = D_q_s(q)

  D(1,1)=(34029*cos(q(4)))/3125 - (1024*cos(q(2) - q(3) + q(4) - q(5)))/3125 - (5971*cos(q(2) -...
          q(3) + q(4)))/3125 - (1024*cos(q(3) - q(2) + q(5)))/3125 - (16*cos(q(2) + q(4)))/5 - (16*cos(q(2)))/5 -...
          (5971*cos(q(2) - q(3)))/3125 + (1024*cos(q(5)))/3125 + 4772949/250000;
  D(1,2)=(34029*cos(q(4)))/3125 - (512*cos(q(2) - q(3) + q(4) - q(5)))/3125 - (5971*cos(q(2) -...
          q(3) + q(4)))/6250 - (512*cos(q(3) - q(2) + q(5)))/3125 - (8*cos(q(2) + q(4)))/5 - (8*cos(q(2)))/5 -...
          (5971*cos(q(2) - q(3)))/6250 + 6658349/500000;
  D(1,3)=(1024*cos(q(5)))/3125 - (512*cos(q(2) - q(3) + q(4) - q(5)))/3125 - (5971*cos(q(2) - q(3) +...
          q(4)))/6250 - (512*cos(q(3) - q(2) + q(5)))/3125 - (5971*cos(q(2) - q(3)))/6250 + 1377549/500000;
  D(1,4)=(34029*cos(q(4)))/6250 - (5971*cos(q(2) - q(3) + q(4)))/6250 - (8*cos(q(2) + q(4)))/5 -...
          (512*cos(q(2) - q(3) + q(4) - q(5)))/3125 + 2204609/312500;
  D(1,5)=(512*cos(q(5)))/3125 - (512*cos(q(3) - q(2) + q(5)))/3125 - (512*cos(q(2) - q(3) + q(4) -...
          q(5)))/3125 + 307009/312500;
  D(2,1)=(34029*cos(q(4)))/3125 - (512*cos(q(2) - q(3) + q(4) - q(5)))/3125 - (5971*cos(q(2) -...
          q(3) + q(4)))/6250 - (512*cos(q(3) - q(2) + q(5)))/3125 - (8*cos(q(2) + q(4)))/5 - (8*cos(q(2)))/5 -...
          (5971*cos(q(2) - q(3)))/6250 + 6658349/500000;
  D(2,2)=(34029*cos(q(4)))/3125 + 6658349/500000;
  D(2,3)=- (5971*cos(q(2) - q(3)))/6250 - (512*cos(q(2) - q(3) + q(4) - q(5)))/3125 - (5971*...
         cos(q(2) - q(3) + q(4)))/6250 - (512*cos(q(3) - q(2) + q(5)))/3125;
  D(2,4)=(34029*cos(q(4)))/6250 + 2204609/312500;
  D(2,5)=- (512*cos(q(2) - q(3) + q(4) - q(5)))/3125 - (512*cos(q(3) - q(2) + q(5)))/3125;
  D(3,1)=(1024*cos(q(5)))/3125 - (512*cos(q(2) - q(3) + q(4) - q(5)))/3125 - (5971*cos(q(2) - q(3) +...
          q(4)))/6250 - (512*cos(q(3) - q(2) + q(5)))/3125 - (5971*cos(q(2) - q(3)))/6250 + 1377549/500000;
  D(3,2)=- (5971*cos(q(2) - q(3)))/6250 - (512*cos(q(2) - q(3) + q(4) - q(5)))/3125 - (5971*...
         cos(q(2) - q(3) + q(4)))/6250 - (512*cos(q(3) - q(2) + q(5)))/3125;
  D(3,3)=(1024*cos(q(5)))/3125 + 1377549/500000;
  D(3,4)=- (512*cos(q(2) - q(3) + q(4) - q(5)))/3125 - (5971*cos(q(2) - q(3) + q(4)))/6250;
  D(3,5)=(512*cos(q(5)))/3125 + 307009/312500;
  D(4,1)=(34029*cos(q(4)))/6250 - (5971*cos(q(2) - q(3) + q(4)))/6250 - (8*cos(q(2) + q(4)))/5 -...
          (512*cos(q(2) - q(3) + q(4) - q(5)))/3125 + 2204609/312500;
  D(4,2)=(34029*cos(q(4)))/6250 + 2204609/312500;
  D(4,3)=- (512*cos(q(2) - q(3) + q(4) - q(5)))/3125 - (5971*cos(q(2) - q(3) + q(4)))/6250;
  D(4,4)=2204609/312500;
  D(4,5)=-(512*cos(q(2) - q(3) + q(4) - q(5)))/3125;
  D(5,1)=(512*cos(q(5)))/3125 - (512*cos(q(3) - q(2) + q(5)))/3125 - (512*cos(q(2) - q(3) + q(4) -...
          q(5)))/3125 + 307009/312500;
  D(5,2)=- (512*cos(q(2) - q(3) + q(4) - q(5)))/3125 - (512*cos(q(3) - q(2) + q(5)))/3125;
  D(5,3)=(512*cos(q(5)))/3125 + 307009/312500;
  D(5,4)=-(512*cos(q(2) - q(3) + q(4) - q(5)))/3125;
  D(5,5)=307009/312500;

 