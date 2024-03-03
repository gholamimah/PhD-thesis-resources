 function seqf = gseq (arraysize)   %% Zigzag function
%  clear all
%  clc
% arraysize = 5
 n = (arraysize + 1)/2;
 arraysize = 2*n-1;
 
 sequence = zeros(2, arraysize^2);
 sequence (1,1) = n;
 sequence (2,1) = n;
 
 dx = +1;
 dy = -1;
 stepx = +1;
 stepy = -1;
 direction = +1;
 counter = 0;
 
 
 for i = 2:arraysize^2
     
     counter = counter +1;
     %%   Direction in the positive sign, movement to forward
     if (direction == +1)
         sequence (1,i) = sequence (1, i-1) +dx;
         sequence (2,i) = sequence (2, i-1);
         
          if (counter == abs (stepx))
              counter = 0;
              direction = -1 * direction;
              dx = -1*dx;
              stepx = -1* stepx;
              
              if stepx >0
                  stepx = stepx +1;
              else 
                  stepx = stepx -1;
              end
          end
     end
     %% 
     
         sequence (1,i) = sequence (1, i-1) ;
         sequence (2,i) = sequence (2, i-1)+ dy;
         
         if (counter == abs (stepy))
              counter = 0;
              direction = -1 * direction;
              dy = -1*dy;
              stepy = -1* stepy;
              
              if stepy >0
                  stepy = stepy +1;
              else 
                  stepy = stepy -1;
              end
         end
 end
         seq = (sequence (1, :) -1) * arraysize +sequence (2,:);
         seqf(1,1:arraysize ^2) = seq;
  end
     