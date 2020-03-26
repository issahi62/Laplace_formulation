% *****************************
%% Finding the capacitance of coaxial rectangle 
%  Developed by: Ibrahim Issah 
%  University of Eastern Finland
%  
%******************************
%% INITIALIZE 
 close all 
 clear all 
 clc
 %****************************
 Nswp = 100;
 Niter = 1000; 
 tolerance = 1e-9; 
 counter = 1; 
 innerwidth = 2;
 innerheight = 2; 
 outerheight =5; 
 outerwidth = 5; 
 
 for sweep = 1:Nswp
 cap = caplap(innerwidth,innerheight,outerwidth,outerheight,sweep,tolerance,Niter);
 M(1, counter) = cap; 
 counter=counter+1; 
 end 
 plot(0:sweep-1, M,'bo', 'Linewidth', 2.5);
 title('convergence values'); 
 xlabel('# of iterations')
 ylabel ('capacitance converging values'); 
 %**************************
 %% CAPACITANCE FORMULATION USING LAPLACE EQUATION DELTA^2U=0;
 %**************************
 
 function cap = caplap(a, b, c, d, n, tol, N)
 %a = width of the inner rectangle 
 %b = height of the inner rectangle 
 %c = width of the outer rectangle 
 %d = height of the outer rectangle 
 %tol = tolerance (1nm)
 % value calculated based on input arguments 
 % relaxation parameter.
 % cap = capacitance per unit length (pF/m)
 % N = Number of iterations
 
 
 rel = 2-c/n; % Relaxation parameter 
 h = .5*c/n;  % Grid size
 na  = round(.5*a/h);  % number of segments on a
 x = linspace(0, 0.5*c, n+1); % Number of segments on d
 m = round(.5*d/h); % Number of segments on d
 mb = round(.5*b/h); % Number of segments on b
 y = linspace(0, .5*d, m+1); % grid points along y-axis
 
 
 % potential array and mask array 
 f = zeros(n+1, m+1); 
 mask = ones( n+1, m+1)*rel; 
 oldcap =0; 
 for i =1:na+1
     for j = 1:mb+1
         mask(i, j) = 0; 
         f(i,j) = 1; 
     end
 end
%   figure(1); imagesc(f)
%   figure(2); imagesc(mask); 
%  
 for iter= 1:N 
     f = seidel(f, mask, n, m); %Perform Guass-siedel iteration 
     cap = guass(n, m, f); % computes the capacitance
     if (abs(cap-oldcap)/cap<tol)
         break
     else
         oldcap= cap; 
     end
 end 
%  str = sprintf('N of iteration = %4i', iter); 
%  fprintf(str);
 end 
 
 
function f = seidel(f, mask, n, m)

 %f is a 2D array with solution after Guass siedel iteration 

 for i = 2:n 
     for j = 2:n 
         f(i,j) = f(i,j)+mask(i,j)*...
             (.25*(f(i-1,j) +f(i+1,j)...
             +f(i,j-1)+f(i,j+1)) -f(i,j));
%          
%          W(i,j)=(.25*(f(i-1,j) +f(i+1,j)...
%              +f(i,j-1)+f(i,j+1)) -f(i,j));
     end
 end
 
 %Symmetr on the left boundary condition i-1->i+1

 i =1;
 for j =2:m 
     f(i,j) = f(i,j)+mask(i,j)*...
             (.25*(f(i+1,j) +f(i+1,j)... 
             +f(i,j-1)+f(i,j+1)) -f(i,j)); %Note the - and +
 end 
 j =1;
 for i =2:m 
     f(i,j) = f(i,j)+mask(i,j)*...
             (.25*(f(i-1,j) +f(i+1,j)... 
             +f(i,j+1)+f(i,j+1)) -f(i,j)); %Note the - and +
 end 
end
 
function cap = guass(n, m, f)

q= 0; 
for i =1:n 
 q= q+(f(i,m)+f(i+1, m))*.5;
end 

for j =1:m
 q= q+(f(n,j)+f(n, j+1))*.5;
end 
 cap = q*4; % 4 quadrants 
 cap = cap*8.854187; %epsilon)*1e12 gives answer in pF/m
end  
 