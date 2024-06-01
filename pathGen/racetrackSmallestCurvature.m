close all;
clear all;

%% racetrack
load('racetrack.mat','t_r'); % load right  boundary from *.mat file
load('racetrack.mat','t_l'); % load left boundary from *.mat file


%% plot left right track boundaries
curveR = cscvn(t_r(:,1:2).');
curveL = cscvn(t_l(:,1:2).');

figure('Name','racetrack','NumberTitle','off','Toolbar','figure','MenuBar','none','OuterPosition',[0 -500 460 1100])
axis equal % eqal axis scaling
axis([-50 70 -50 450]) % plot height and width
fnplt(curveR, 'b', 1); 
hold on;
fnplt(curveL, 'b', 1);

%% solver
t_r = t_r(1:1950,1:2).';
t_l = t_l(1:1950,1:2).';

curveR = cscvn(t_r);
curveL = cscvn(t_l);

curveLxCoefs = curveL.coefs(1:2:end-1,:);
curveLyCoefs = curveL.coefs(2:2:end,:);
curveRxCoefs = curveR.coefs(1:2:end-1,:);
curveRyCoefs = curveR.coefs(2:2:end,:);

xL = curveLxCoefs(:,4);
yL = curveLyCoefs(:,4);
xR = curveRxCoefs(:,4);
yR = curveRyCoefs(:,4);

xC = xR + 0.5.*(xL - xR);
yC = yR + 0.5.*(yL - yR);
distance = sqrt(power(xL - xR,2) + power(yL - yR,2));

%%
steps = 2000;
centerLine = [xC yC];
stepLengths = sqrt(sum(diff(centerLine,[],1).^2,2));
stepLengths = [0; stepLengths]; % add the starting point
cumulativeLength = cumsum(stepLengths);
stepLocations = linspace(0,cumulativeLength(end), steps);
centerLine = interp1(cumulativeLength, centerLine, stepLocations);
distance = interp1(cumulativeLength, distance, stepLocations)';
xC = centerLine(:,1);
yC = centerLine(:,2);

% normal direction for each vertex
dx = gradient(xC);
dy = gradient(yC);
dL = hypot(dx,dy);

% Smoothing
%distance = atan(distance-5) + 5;

% Reconstruction of track boundaries
xoff = -0.5*distance.*dy./dL;
yoff = 0.5*distance.*dx./dL;
xL = xC + xoff;      % get inner offset curve
yL = yC + yoff;
xR = xC - xoff;      % get outer offset curve
yR = yC - yoff;

xC = [xC; xC(1,1)];
xL = [xL; xL(1,1)];
yL = [yL; yL(1,1)];
xR = [xR; xR(1,1)];
yR = [yR; yR(1,1)];

% form delta Matrices
delx = xL - xR;
dely = yL - yR;

n = numel(delx);

% preallocation
H = zeros(n);
B = zeros(size(delx)).';

% formation of H matrix (nxn)
for i=2:n-1
    
    % first row
    H(i-1,i-1) = H(i-1,i-1) + delx(i-1)^2         + dely(i-1)^2;
    H(i-1,i)   = H(i-1,i)   - 2*delx(i-1)*delx(i) - 2*dely(i-1)*dely(i);
    H(i-1,i+1) = H(i-1,i+1) + delx(i-1)*delx(i+1) + dely(i-1)*dely(i+1);
    
    %second row
    H(i,i-1)   = H(i,i-1)   - 2*delx(i-1)*delx(i) - 2*dely(i-1)*dely(i);
    H(i,i)     = H(i,i )    + 4*delx(i)^2         + 4*dely(i)^2;
    H(i,i+1)   = H(i,i+1)   - 2*delx(i)*delx(i+1) - 2*dely(i)*dely(i+1);
    
    % third row
    H(i+1,i-1) = H(i+1,i-1) + delx(i-1)*delx(i+1) + dely(i-1)*dely(i+1);
    H(i+1,i)   = H(i+1,i)   - 2*delx(i)*delx(i+1) - 2*dely(i)*dely(i+1);
    H(i+1,i+1) = H(i+1,i+1) + delx(i+1)^2         + dely(i+1)^2;
    
end

% formation of B matrix (1xn)
for i=2:n-1
    
    B(1,i-1) = B(1,i-1) + 2*(xR(i+1)+xR(i-1)-2*xR(i))*delx(i-1) + 2*(yR(i+1)+yR(i-1)-2*yR(i))*dely(i-1);
    B(1,i)   = B(1,i)   - 4*(xR(i+1)+xR(i-1)-2*xR(i))*delx(i)   - 4*(yR(i+1)+yR(i-1)-2*yR(i))*dely(i);
    B(1,i+1) = B(1,i+1) + 2*(xR(i+1)+xR(i-1)-2*xR(i))*delx(i+1) + 2*(yR(i+1)+yR(i-1)-2*yR(i))*dely(i+1);
    
end

% define constraints
lb = zeros(n,1) + 0.25;
ub = ones(size(lb)) - 0.25;

% if start and end points are the same
Aeq      =   zeros(2,n);
Aeq(1,1)   =   1;
Aeq(1,end) =   -1;
Aeq(2,1) = 1;
beq =   [0; 0.5];
    
%% Solver

options = optimoptions('quadprog','Display','iter');
[resMCP,fval,exitflag,output] = quadprog(2*H,B',[],[],Aeq,beq,lb,ub,[],options);

% co-ordinates for the resultant curve
xresMCP = zeros(size(xC));
yresMCP = zeros(size(xC));

for i = 1:numel(xC)
    xresMCP(i) = xR(i)+resMCP(i)*delx(i);
    yresMCP(i) = yR(i)+resMCP(i)*dely(i);
end

plot(xresMCP,yresMCP,'color','r')
axis equal
tMCP = [xresMCP,yresMCP];

save("referenceMCP.mat", "tMCP")