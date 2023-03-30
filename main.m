%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Title: 
% Desc: 
% Author: Jc Chen
% Modified: 2023/03/29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clc; clear;
addpath EchoGenerate/ RDA/ AnalyseTool/;

% 场景参数设置
c = 3e8;
azm_len = 260;
rng_len = 520;
rng_start = 10000;

% 分辨率要求
rho_r = 1;
rho_a = 1;
% 过采样要求
os_r = 1.2;
os_a = 1.2;

% 雷达参数
V = 100;
La = 2 * rho_a;

theta_rc = 0;
theta_rc = 0.00*pi;
% 脉冲相关
f0 = 2e9;
lamda = c/f0;
theta_bw = 0.886*lamda/La;
Tp = 5e-6;
B = 0.886*c/2/rho_r;
Kr = B/Tp;

f_dop = 2*V*cos(theta_rc)/lamda*theta_bw;       % 多普勒带宽
PRF = os_a * f_dop;
fs = os_r * B;
fnc = 2*V*sin(theta_rc)/lamda;


% ****************** 回波生成 **************** %

%tg_pos = [260 0; 261 1];
tg_pos = [260 0];           % 设置点目标xy位置
data = getSimulateEcho(azm_len, rng_len, rng_start, f0, fs, PRF, V, Kr, Tp, theta_rc, theta_bw, tg_pos);

% ****************** 算法成像 **************** %
img = RDA(data, Kr, f0, fs, PRF, V, fnc, rng_start, theta_bw=theta_bw, Tp=Tp);

% ****************** 结果 **************** %
imagesc(abs(img));
axis xy
inspectView(img, 20, [64 40]);


