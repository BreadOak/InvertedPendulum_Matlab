clc
clear all

global M
global m
global L
global g
global y

Start_time = 0;
End_time = 20;
Sampling_Frequency = 30;
st = 1/Sampling_Frequency;
tspan = [0,st];

l = 1 : 1 : (End_time-Start_time)/st;

g = -9.81; % gravity(ms/s^2)
L = 2;     % length of pendulum(m)
m = 1;     % mass of pendulum(kg)
M = 3;     % mass of cart(kg)

y = 0.5;   % hight of cart's center(m)
b = 1;     % pendulum up (b = 1)

A = [0,          1,                0, 0;
     0,          0,          b*m*g/M, 0;
     0,          0,                0, 1;
     0,          0, -b*(M+m)*g/(M*L), 0];

B = [      0;
         1/M;
           0;
      b/(M*L)];

C = [1,0,0,0;
     0,1,0,0;
     0,0,1,0;
     0,0,0,1];
 
D = [0;0;0;0];

Q = [1,0,0,0;
     0,1,0,0;
     0,0,1,0;
     0,0,0,1];

R = 0.0001;
u = 0;

%% Initial state
x_state1_ini = 0;        %(m)
x_state2_ini = 0;        %(m/s)
x_state3_ini = pi + pi/9;%(rad)
x_state4_ini = 0;        %(rad/s)
x_state_ini = [x_state1_ini;x_state2_ini;x_state3_ini;x_state4_ini];

%% Desired state
x_state1_des = 10; %(m)
x_state2_des = 0;  %(m/s)
x_state3_des = pi; %(rad)
x_state4_des = 0;  %(rad/s)
x_state_des = [x_state1_des;x_state2_des;x_state3_des;x_state4_des];

%% Reset initial data
x_data = [];
dx_data = [];
th_data = [];
dth_data = [];

x_err_data = [];
dx_err_data = [];
th_err_data = [];
dth_err_data = [];

d_hat_data = [];
d_data = [];

%% Calculate LQR gain
[K,P] = lqr(A,B,Q,R);
LQR_gain = K;

%% Control loop
for i = l    
    x_error = x_state_ini - x_state_des;
    
    u = -LQR_gain*x_error;
    
    % plant
    [t,S_list] = ode45(@(t,x_state_ini) InvertedPendulum(t,x_state_ini,u), tspan, x_state_ini);
    state = S_list(end,(1:4)).';
    
    % sate update
    x_state_ini = state;
    
    % save data
    x_data(end+1) = state(1);
    dx_data(end+1) = state(2);
    th_data(end+1) = state(3);
    dth_data(end+1) = state(4);
    
    x_err_data(end+1) = x_error(1);
    dx_err_data(end+1) = x_error(2);
    th_err_data(end+1) = x_error(3);
    dth_err_data(end+1) = x_error(4);
end

%% Simulation
x = x_data;
th = th_data;

for point = l
    draw_robot(x(point),th(point));
    %pause(0.01)
end

%% Error plot
f1 = figure;
subplot(2,1,1)
plot(l*st,x_err_data)
xlabel('time')
ylabel('X error')
xlim([0, End_time])
ylim([-15, 15])
grid on

subplot(2,1,2)
plot(l*st,th_err_data)
xlabel('time')
ylabel('Theta error')
xlim([0, End_time])
ylim([-pi, pi])
grid on

f2 = figure;
subplot(2,1,1)
plot(l*st,dx_err_data)
xlabel('time')
ylabel('dX error')
xlim([0, End_time])
ylim([-20, 20])
grid on

subplot(2,1,2)
plot(l*st,dth_err_data)
xlabel('time')
ylabel('dTheta error')
xlim([0, End_time])
ylim([-3*pi, 3*pi])
grid on

%% Inverted Pendulum dynamics
function dxdt = InvertedPendulum(t,x,u)
    global M
    global m
    global g
    global L
    D = m*L^2*(M + m*(1-cos(x(3))^2));
    dx_dt1 = x(2);
    dx_dt2 = (1/D)*(-m^2*L^2*g*cos(x(3))*sin(x(3)) + m*L^2*(m*L*x(4)^2*sin(x(3)))) + m*L^2*(1/D)*u;
    dx_dt3 = x(4);
    dx_dt4 = (1/D)*((M+m)*m*g*L*sin(x(3)) - m*L*cos(x(3))*(m*L*x(4)^2*sin(x(3)))) - m*L*cos(x(3))*(1/D)*u;
    dxdt = [dx_dt1;dx_dt2;dx_dt3;dx_dt4];
end

%% Draw robot
function draw_robot(x,theta)
    global L
    global y
    plot(0,-1,'k.')
    hold on;
    grid on;
    line([x,x+L*sin(theta)],[y,-L*cos(theta)+y])
    drawCircle(x+L*sin(theta),-L*cos(theta)+y,0.2)
    rectangle('Position',[x-1 0 2 1])
    axis([-15,15,-1,5])
    daspect([1,1,1])
    drawnow;
    hold off;
end

%% Draw Circle
function drawCircle(x,y,r)
    th = linspace(0,2*pi,10);
    plot(x+r*cos(th),y+r*sin(th),'r')
end