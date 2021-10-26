%{
MIT License

Copyright (c) 2021 ml0pe5

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

%}


% Dentate gyrus model, modified from [1]
% [1] López-Cuevas et al., NeuroImage 113, 2015

clc; clear; close all;

% Parameters (see Table S37 for their meaning)
A  = 5;
B  = 10;
G  = 10;

a = 100; 
b = 50; 
g = 500;  

x1 = 10; x2 = 0.5; x3 = 1;
C1 = 16*x1;   % MC to GC
C2 = 1*x1;    % BC to GC
C3 = 8*x1;    % HC to GC
C4 = 400*x2;  % GC to MC
C5 = 20*x3;   % GC to BC
C6 = 80*x3;   % MC to BC
C7 = 100;     % GC to HC

MC_absent = 0; % set to 1 to remove MC
if MC_absent
    C1 = 0; % to remove MC to GC
    C6 = 0; % to remove MC to BC
end

p = 0.2; % sigma = p/sqrt(dt)

T = 20000; 
dt = 0.001;         

% Stimulus:
s = linspace(0,1.5,T); % increasing stimulus 
% s = zeros(T,1); % no stimulus 
% s = 1.1*ones(T,1); % DC stimulus


% Initial conditions:
v_gc = zeros(T,1);  x_gc = 0;
v_mc = zeros(T,1);  x_mc = 0;
v_bc = zeros(T,1);  x_bc = 0;
v_hc = zeros(T,1);  x_hc = 0;
time = zeros(T,1);

for t = 2:T    
    dv_gc = x_gc; 
    dv_mc = x_mc; 
    dv_bc = x_bc; 
    dv_hc = x_hc; 
     
    dx_gc = A*a*(p*randn+s(t-1)+S_fun(C1*v_mc(t-1)-C2*v_bc(t-1)-...
        C3*v_hc(t-1)))-2*a*x_gc-a^2*v_gc(t-1);
    dx_mc = A*a*S_fun(C4*v_gc(t-1))-2*a*x_mc-a^2*v_mc(t-1);
    dx_bc = G*g*S_fun(C5*v_gc(t-1)+C6*v_mc(t-1))-2*g*x_bc-g^2*v_bc(t-1);
    dx_hc = B*b*(p*randn+S_fun(C7*v_gc(t-1)))-2*b*x_hc-b^2*v_hc(t-1);
    
    v_gc(t) = v_gc(t-1) + dt*dv_gc;
    v_mc(t) = v_mc(t-1) + dt*dv_mc;
    v_bc(t) = v_bc(t-1) + dt*dv_bc;
    v_hc(t) = v_hc(t-1) + dt*dv_hc;
    
    x_gc = x_gc + dt*dx_gc;
    x_mc = x_mc + dt*dx_mc;
    x_bc = x_bc + dt*dx_bc;
    x_hc = x_hc + dt*dx_hc;        
    
    time(t) = time(t-1)+dt;
end

% Plot results: 
hfig = figure;
plot(time,v_gc,'k','linewidth',2)
ylabel('v_{gc}')
xlabel('time')



% Sigmoid function S:
function S=S_fun(v)
% Parameters (see Table S37 for their meaning)
v_0=6;
e_0=2.5;
r=0.56;

S=2*e_0/(1+exp(r*(v_0-v)));
end

