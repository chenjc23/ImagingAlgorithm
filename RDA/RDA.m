function RDA()
global azm_len, rng_len, rng_start, V, PRF, theta_rc, La, f0, Tp, B, fs
global theta_bw
global sra
global data

lamda = c / f0;
fnc = 2*V*sin(theta_rc)/lamda;



Ls = theta_bw * (rng_start+rng_len) / cos(theta_rc);        % 合成孔径长度
Ta = Ls / V;                                                % 合成孔径时间

global ta, tr;
%tr = -Tp/2 + rng_start*2/c : 1/fs : (rng_start+rng_len)*2/c + Tp/2;        % 距离向时间序列

Saf = fft(data, [], 2);                   % 距离频域
Sdf = fft(Saf, [], 1);                    % 二维频域

[Na, Nr] = size(data);

fr = 0:fs/Nr:fs-1/Nr;                     % 距离向频率序列
fa = circshift(0: PRF/Na: PRF - PRF/Na, fix(fnc*Na/PRF)).';         % 按多普勒中心进行循环移位构建多普勒序列


D = sqrt(1-lamda^2*fa.^2/4/V^2);           % 徙动因子
R_ref = rng_start + rng_len/2;                                          % 测绘带中心参考距离
Ksrc = 2*V^2*f0^3*D^3 ./ (c*R_ref*fa.^2);      % 二次调频率
Hm = exp(1j*pi*(1/Kr-1./Ksrc) * fr.^2);                               % 一次二次距离压缩滤波器

Srd = ifft(Sdf .* Hm, [], 2);             % 距离压缩后变换到rd域

% *********** RCMC ************ %
Rn = c/2 * tr;            % 距离序列
del_R = lamda^2* (fa.^2) * Rn/8/V^2;      % RCM
core_len = 8;
Srd_RCMC = zeros(Na, Nr);
for i = 1:Na
  for j = 1:Nr
    ideal_p = j + 2*del_R(i,j)/c / (1/fs);          % 需要插值的浮点位置
    quant_p = fix(ideal_p);                         % 需要插值的量化位置
    if (quant_p < core_len)


  end
end




Haz = exp(1j*4*pi/lamda * D.' * Rn);      % 方位压缩滤波器

















end