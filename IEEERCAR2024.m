%% Research code by Agus Hasan
%% This code is used for numerical simulation of flexible joint robot system using WyNDA algorithm

clear;
clc;

%% time horizon
tf  = 20;           % define the simulation horizon
dt  = 0.001;        % define the time-step (the sampling frequency)
t   = dt:dt:tf;     % define the time array

%% number of variables and coefficients
n = 4;              % define the number of measured state
r = 60;             % define the number of basis function

%% noise
R = 0.001;          % define the standard deviation of measurement noise

%% state initialization
x        = zeros(n,1);      % define the initial value of the robot
xbar     = x;               % define the initial value of the state estimates
y        = x;               % define the initial value of the measurement
thetabar = zeros(r,1);      % define the initial value of the parameter estimate
 
%% true parameters
Ks      = 1.61;
Jh      = 0.0021;
M       = 0.403;
g       = -0.98;
h       = 0.06;
Km      = 0.00767;
Kg      = 70;
Jl      = 0.0059;
Rm      = 2.6;
% augmented parameters
alpha   = Ks/Jh;
beta    = (Km^2*Kg^2)/(Rm*Jh);
gamma   = Ks/Jl;
delta   = M*g*h/Jl;
omega   = (Km*Kg)/(Rm*Jh);

%% initial control inputs
u     = 0;

%% for plotting
uArray          = [];
xArray          = [];
xbarArray       = [];
yArray          = [];
thetabarArray   = [];

%% Initialization for the adaptive observer
lambdav = 0.999;
lambdat = 0.999;
Rx      = 0.1*eye(n);
Rt      = 1*eye(n);
Px      = 1000*eye(n);
Pt      = 10*eye(r);
Gamma   = 0*ones(n,r);

%% simulation
for i=1:(tf/dt)

    u = 4*cos(2*i*dt);

    uArray         = [uArray u];
    xArray         = [xArray x];
    xbarArray      = [xbarArray xbar];    
    yArray         = [yArray y];
    thetabarArray  = [thetabarArray thetabar]; 

    %% Simulate the system using Runge-Kutta
    % the objective of this simulation is to generate data from the robot
    % system
    k1 = x(3);
    l1 = x(4);
    m1 = alpha*x(2)-beta*x(3)+omega*u;
    n1 = -(alpha+gamma)*x(2)+delta*sin(x(1)+x(2))+beta*x(3)-omega*u;
    k2 = x(3)+0.5*dt*m1;
    l2 = x(4)+0.5*dt*n1;
    m2 = alpha*(x(2)+0.5*dt*l1)-beta*(x(3)+0.5*dt*m1)+omega*u;
    n2 = -(alpha+gamma)*(x(2)+0.5*dt*l1)+delta*sin((x(1)+0.5*dt*k1)+(x(2)+0.5*dt*l1))+beta*(x(3)+0.5*dt*m1)-omega*u;
    k3 = x(3)+0.5*dt*m2;
    l3 = x(4)+0.5*dt*n2;
    m3 = alpha*(x(2)+0.5*dt*l2)-beta*(x(3)+0.5*dt*m2)+omega*u;
    n3 = -(alpha+gamma)*(x(2)+0.5*dt*l2)+delta*sin((x(1)+0.5*dt*k2)+(x(2)+0.5*dt*l2))+beta*(x(3)+0.5*dt*m2)-omega*u;
    k4 = x(3)+dt*m3;
    l4 = x(4)+dt*n3;
    m4 = alpha*(x(2)+dt*l3)-beta*(x(3)+dt*m3)+omega*u;
    n4 = -(alpha+gamma)*(x(2)+dt*l3)+delta*sin((x(1)+dt*k3)+(x(2)+dt*l3))+beta*(x(3)+dt*m3)-omega*u;
    x(1) = x(1) + (dt/6)*(k1+2*k2+2*k3+k4);
    x(2) = x(2) + (dt/6)*(l1+2*l2+2*l3+l4);
    x(3) = x(3) + (dt/6)*(m1+2*m2+2*m3+m4);
    x(4) = x(4) + (dt/6)*(n1+2*n2+2*n3+n4);

    y = x+dt*R^2*randn(n,1);
    
    %% here we define the approximation function, which consists of as set of basis function

    Phi = [y(1) y(2) y(3) y(4) y(1)^2 y(2)^2 y(3)^2 y(4)^2 sin(y(1)+y(2)) sin(y(1)+y(3)) sin(y(1)+y(4)) sin(y(2)+y(3)) sin(y(2)+y(4)) sin(y(3)+y(4)) u zeros(45,1)';
           zeros(15,1)' y(1) y(2) y(3) y(4) y(1)^2 y(2)^2 y(3)^2 y(4)^2 sin(y(1)+y(2)) sin(y(1)+y(3)) sin(y(1)+y(4)) sin(y(2)+y(3)) sin(y(2)+y(4)) sin(y(3)+y(4)) u zeros(30,1)';
           zeros(30,1)' y(1) y(2) y(3) y(4) y(1)^2 y(2)^2 y(3)^2 y(4)^2 sin(y(1)+y(2)) sin(y(1)+y(3)) sin(y(1)+y(4)) sin(y(2)+y(3)) sin(y(2)+y(4)) sin(y(3)+y(4)) u zeros(15,1)';
           zeros(45,1)' y(1) y(2) y(3) y(4) y(1)^2 y(2)^2 y(3)^2 y(4)^2 sin(y(1)+y(2)) sin(y(1)+y(3)) sin(y(1)+y(4)) sin(y(2)+y(3)) sin(y(2)+y(4)) sin(y(3)+y(4)) u];
    % if the structure of the governing equation is known, we can use the
    % following approximation function. The problem becomes parameter
    % identification problem
    % do not forget to change the number of basis function r
%    Phi = [y(3) 0 0 0 0 0 0 0 0;
%           0 y(4) 0 0 0 0 0 0 0;
%           0 0 y(2) y(3) u 0 0 0 0;
%           0 0 0 0 0 y(2) sin(y(1)+y(2)) y(3) u];    
    
    %% Estimation using adaptive observer
    % these lines contain the codes for adaptive observer
    Kx = Px*inv(Px+Rx);
    Kt = Pt*Gamma'*inv(Gamma*Pt*Gamma'+Rt);
    Gamma = (eye(n)-Kx)*Gamma;

    xbar = xbar+(Kx+Gamma*Kt)*(y-xbar);
    thetabar = thetabar-Kt*(y-xbar);

    xbar = xbar+Phi*thetabar;

    thetabar = thetabar;
    Px = (1/lambdav)*eye(n)*(eye(n)-Kx)*Px*eye(n);
    Pt = (1/lambdat)*(eye(r)-Kt*Gamma)*Pt;
    Gamma = eye(n)*Gamma-Phi;

end

figure(1)
subplot(2,2,1)
plot(t,yArray(1,:),'b-','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',56)
legend('measurement')
ylabel('p')
grid on;
grid minor;
subplot(2,2,2)
plot(t,yArray(2,:),'b-','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',56)
ylim([-0.1 0.1])
ylabel('q')
grid on;
grid minor;
subplot(2,2,3)
plot(t,yArray(3,:),'b-','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',56)
ylabel('r')
xlabel('t (s)')
grid on;
grid minor;
subplot(2,2,4)
plot(t,yArray(4,:),'b-','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',56)
ylim([-0.4 0.4])
ylabel('s')
xlabel('t (s)')
grid on;
grid minor;

figure(2)
subplot(2,2,1)
plot(t,yArray(1,:),'b-','LineWidth',10);
hold on;
plot(t,xbarArray(1,:),'g:','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',56)
legend('measured','estimated')
ylabel('p')
grid on;
grid minor;
subplot(2,2,2)
plot(t,yArray(2,:),'b-','LineWidth',10);
hold on;
plot(t,xbarArray(2,:),'g:','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',56)
ylim([-0.1 0.1])
ylabel('q')
grid on;
grid minor;
subplot(2,2,3)
plot(t,yArray(3,:),'b-','LineWidth',10);
hold on;
plot(t,xbarArray(3,:),'g:','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',56)
ylabel('r')
xlabel('t (s)')
grid on;
grid minor;
subplot(2,2,4)
plot(t,yArray(4,:),'b-','LineWidth',10);
hold on;
plot(t,xbarArray(4,:),'g:','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',56)
ylim([-0.4 0.4])
ylabel('s')
xlabel('t (s)')
grid on;
grid minor;

figure(3)
stem([1:1:60],log(abs(thetabarArray(:,tf/dt)/dt)),'LineWidth',4,'LineStyle','-.',...
     'MarkerFaceColor','red',...
     'MarkerEdgeColor','green','MarkerSize',12);
set(gca,'color','white','LineWidth',3,'FontSize',72)
grid on;
grid minor;
ylabel('log(\theta_i)','FontSize',72)
xlabel('i','FontSize',72)

figure(4)
subplot(3,2,1)
plot(t,uArray,'-k','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
legend('control input')
grid on;
grid minor;
ylabel('u','FontSize',72)
subplot(3,2,2)
plot(t,alpha*ones(1,length(t)),'b-','LineWidth',10);
hold on;
plot(t,thetabarArray(32,:)/dt,'g:','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
legend('true parameter','estimated parameter')
ylim([0 1000])
grid on;
grid minor;
ylabel('\alpha','FontSize',72)
subplot(3,2,3)
plot(t,beta*ones(1,length(t)),'b-','LineWidth',10);
hold on;
plot(t,-thetabarArray(33,:)/dt,'g:','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
ylim([-200 200])
grid on;
grid minor;
ylabel('\beta','FontSize',72)
subplot(3,2,4)
plot(t,gamma*ones(1,length(t)),'b-','LineWidth',10);
hold on;
plot(t,-thetabarArray(47,:)/dt-thetabarArray(32,:)/dt,'g:','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
ylim([100 400])
grid on;
grid minor;
ylabel('\gamma','FontSize',72)
subplot(3,2,5)
plot(t,delta*ones(1,length(t)),'b-','LineWidth',10);
hold on;
plot(t,thetabarArray(54,:)/dt,'g:','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
ylim([-400 400])
grid on;
grid minor;
ylabel('\delta','FontSize',72)
xlabel('t (s)')
subplot(3,2,6)
plot(t,omega*ones(1,length(t)),'b-','LineWidth',10);
hold on;
plot(t,thetabarArray(45,:)/dt,'g:','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
ylim([0 200])
grid on;
grid minor;
ylabel('\omega','FontSize',72)
xlabel('t (s)')