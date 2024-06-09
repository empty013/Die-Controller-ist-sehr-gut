function [Y, U_arr, debug] = ode1_s(odefun,sspan,y0,varargin)
%ODE1  Solve differential equations with a non-adaptive method of order 1.
%   Y = ODE1(ODEFUN,TSPAN,Y0) with TSPAN = [T1, T2, T3, ... TN] integrates 
%   the system of differential equations y' = f(t,y) by stepping from T0 to 
%   T1 to TN. Function ODEFUN(T,Y) must return f(t,y) in a column vector.
%   The vector Y0 is the initial conditions at T0. Each row in the solution 
%   array Y corresponds to a time specified in TSPAN.
%
%   Y = ODE1(ODEFUN,TSPAN,Y0,P1,P2...) passes the additional parameters 
%   P1,P2... to the derivative function as ODEFUN(T,Y,P1,P2...). 
%
%   This is a non-adaptive solver. The step sequence is determined by TSPAN.
%   The solver implements the forward Euler method of order 1.   
%
%   Example 
%         tspan = 0:0.1:20;
%         y = ode1(@vdp1,tspan,[2 0]);  
%         plot(tspan,y(:,1));
%     solves the system y' = vdp1(t,y) with a constant step size of 0.1, 
%     and plots the first component of the solution.   
%

if ~isnumeric(sspan)
  error('TSPAN should be a vector of integration steps.');
end

if ~isnumeric(y0)
  error('Y0 should be a vector of initial conditions.');
end

h = diff(sspan);
if any(sign(h(1))*h <= 0)
  error('Entries of TSPAN are not in order.') 
end  

% try
%   f0 = feval(odefun,sspan(1),y0,varargin{:});
% catch
%   msg = ['Unable to evaluate the ODEFUN at t0,y0. ',lasterr];
%   error(msg);  
% end  

% y0 = y0(:);   % Make a column vector.
% if ~isequal(size(y0),size(f0))
%   error('Inconsistent sizes of Y0 and f(t0,y0).');
% end  

neq = length(y0);
N = length(sspan);
Y = zeros(neq,N);
U_arr = zeros(2, N-1);
debug = zeros(4, N);

% q = [e_y; e_psi; v; beta; omega; t];
[~, psi_ref, x_ref, y_ref] = referencePath(0);
Y(:,1) = y0;
% x_ref = -2.5;
% y_ref = 0;
% psi_ref = pi/2;
debug(1, 1) = x_ref - y0(1) * sin(psi_ref);
debug(2, 1) = y_ref + y0(1) * cos(psi_ref);
debug(3, 1) = x_ref;
debug(4, 1) = y_ref;
tic
for i = 1:N-1 
  % i
  if mod(i, 10) == 0
      i
  end
  [dY_i, U_i] = feval(odefun,sspan(i),Y(:,i),varargin{:});
  Y(:,i+1) = Y(:,i) + h(i)*dY_i;
  [~, psi_ref, x_ref, y_ref] = referencePath(sspan(i+1));
  debug(1, i+1) = x_ref - Y(1, i+1) * sin(psi_ref);
  debug(3, i+1) = x_ref;
  % y_ref = sspan(i+1);
  debug(2, i+1) = y_ref + Y(1, i+1) * cos(psi_ref);
  debug(4, i+1) = y_ref;
  U_arr(:,i) = U_i;
  % debug(:,i) = debug_i;
end
toc
Y = Y.';
U_arr = U_arr.';
debug = debug.';
