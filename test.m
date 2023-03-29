%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Title: 
% Desc: 
% Author: Jc Chen
% Modified: 2023/03/29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clc; clear;
addpath EchoGenerate/ RDA/ AnalyseTool/ data/;

load CDdata1.mat data;
load CD_run_params.mat ...
  Fr Kr PRF R0 Tr c f0 Nrg_cells Nrg_lines_blk first_rg_cell first_rg_line;         % 成像处理相关的参数

Nr = Nrg_cells;                 % 距离门数目
Na = Nrg_lines_blk;             % 距离线数目

lamda = c/f0;                  % 信号波长
Kr = -Kr;                       % 发射脉冲负扫频
V = 7062;                      % 直线几何约束下有效雷达速率（忽略空变，即为参考速度）
fnc = -6900;                   % 多普勒中心频率
Tp = Tr;
fs = Fr;

rng_start = R0 - 1/fs*c/2 * fix(Nr/2);

% ****************** 算法成像 **************** %
img = RDA(data, Kr, f0, Tp, fs, PRF, V, rng_start, fnc=fnc, type="scence");

% ****************** 结果 **************** %
s_enhance = 20*log10(abs(img)/max(max(abs(img)))+eps);
imagesc(s_enhance, [-65, 0]);
colormap("gray");
axis xy



