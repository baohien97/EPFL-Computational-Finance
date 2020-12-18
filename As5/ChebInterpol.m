function [ IN ] = ChebInterpol(f, x, n, a, b)
    % Chebyshev interpolation of Function f for parameter x of degree n3
    % 
    % % Input:
    % f: Function to interpolate (function handle)
    % x: Points to be interpolated
    % n: Number of Chebyhsev Points (can be multidim)
    % [a,b]: Interval of Parameter (optional)
   %  Output: IN: interpolated prices1011
   
  %  Set interval if necessary (Chebyshev points will be set on this interval)
  if ~exist('a','var')
      a = min(x);
  end
  if ~exist('b','var')
      b = max(x);
  end
  % Initiate Interpolation output
  IN = zeros(length(x),length(n));
  % linear tranformation from larger intervals to [-1,1]
  x = (x-a)./(b-a)*2-1;
  for j = 1:length(n) 
      %to calculate Interpolation for different Ns at once
     %  Set Chebyshev points
     k = 0:n(j);
     pp = cos(pi*(k)/(n(j)));
    %  Linear tranformation for larger intervals
    pp = (b-a)/2*(pp+1)+a;
    % Initiate weights c
    c = zeros(1,n(j)+1);
    try
        PriceChebpts = f(pp);
    catch
        PriceChebpts = zeros(size(pp));
        for i = 1:size(pp,2)
            PriceChebpts(i) = f(pp(i));
        end    
    end
    for m = 0:n(j)
        c(m+1) = ...
            2.^(m>0)/n(j)*(PriceChebpts*diag([0.5;ones(n(j)-1,1);0.5])*...
            (cos(m*pi*(k)/(n(j))))'); %diag to take first an last summand 
            %...times 0.5
    end
    % Interpolated price
    for i = 1:length(x)
        T = cos(k*acos(x(i)));
        IN(i,j) = dot(c,T);
    end
  end
end
  
