
F = ...
   [0 0 0 0;...
    1 1 0 0;...
    0 0 0 0;...
    2 2 0 0];

prematrix = ...
    [+1 +0 +0 +0;...
     +0 +0 +1 +0;...
     -3 +3 -2 -1;...
     +2 -2 +1 +1];
 
postmatrix = ...
    [+1 +0 -3 +2;...
     +0 +0 +3 -2;...
     +0 +1 -2 +1;...
     +0 +0 -1 +1];
 
 alpha = prematrix * F * postmatrix
 
 prematrix2 = ...
    [+1 +0 +0 +0;...
     +1 +1 +1 +1;...
     +0 +1 +0 +0;...
     +0 +1 +2 +3];
 
 postmatrix2 = ...
    [+1 +1 +0 +0;...
     +0 +1 +1 +1;...
     +0 +1 +0 +2;...
     +0 +1 +0 +3];
 
 F_recovered = prematrix2 * alpha * postmatrix2;
 
 F == F_recovered;
 
 x_test = [0:0.01:1];
 y = 0;
 z = zeros(length(x_test),1);
 
 for i = 1:length(x_test)
     x = x_test(i)/1;
     x_vector = [1 x x^2 x^3];
     y_vector = [1; y; y^2; y^3];
     
     z(i) = x_vector * alpha * y_vector;
 end
 
 figure; plot(x_test, z)
 
 hold on;
 plot(x_test, x_test.^2, 'r')
 
 
%  x = 0;
%  y_test = [0:0.01:1];
%  z = zeros(length(y_test),1);
%  
%  for i = 1:length(y_test)
%      y = y_test(i);
%      x_vector = [1 x x^2 x^3];
%      y_vector = [1; y; y^2; y^3];
%      
%      z(i) = x_vector * alpha * y_vector;
%  end
%  
%  figure; plot(y_test, z)
     
     
     
     
     
     