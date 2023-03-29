function RDA(data, fs, PRF, Kr, f0, Tp, theta_rc, theta_bw, V, rng_start)
c = 3e8;
[na, nr] = size(data);
rng_len = 1/fs*c/2 * (nr-1);



lamda = c / f0;
fnc = 2*V*sin(theta_rc)/lamda;


Ls = theta_bw * (rng_start+rng_len) / cos(theta_rc);        % 合成孔径长度
Ta = Ls / V / cos(theta_rc);                                % 合成孔径时间


tr = rng_start*2/c : 1/fs : (rng_start+rng_len)*2/c + 1.2*Tp;        % 距离向时间序列

Na = fix(Ta*PRF) + na;
Nr = length(tr);


%sra = zeros(Na, Nr);

Saf = fft(data, Nr, 2);                   % 距离频域
Sdf = fft(Saf, Na, 1);                    % 二维频域


                  
fr = fs/Nr * ((0:Nr-1)-fix(Nr/2));            % 距离向频率序列
fr = circshift(fr, -fix(Nr/2));            % 零频移到两端

fa = PRF/Na * ((0:Na-1)-fix(Na/2)).';         
fa = circshift(fa, -fix(Na/2)) + fnc;      % 按多普勒中心进行循环移位构建多普勒序列



R_ref = rng_start + rng_len/2;                 % 测绘带中心参考距离
D = sqrt(1-lamda^2*fa.^2/4/V^2);               % 徙动因子
Ksrc = 2*V^2*f0^3*D.^3 ./ (c*R_ref*fa.^2);      % 二次调频率
Hr = exp(1j*pi*(1/Kr-1./Ksrc) * fr.^2);        % 一次二次距离压缩滤波器
Srd = ifft(Sdf .* Hr, [], 2);                  % 距离压缩后变换到rd域

Srd = Srd(:,1:nr);                             % 距离向舍弃弃置区
Nr = nr;
tr = tr(1:nr);

shiyu = ifft(Srd, [], 1);
temp = abs(shiyu(fix(na/2), :));
temp = 20*log10(temp/max(temp));

figure
plot(temp)




figure
contour(abs(ifft(Srd, [], 1)));
figure
imagesc(abs(Srd))

% *********** RCMC ************ %
Rn = c/2 * tr;            % 距离序列
del_R = lamda^2* (fa.^2) * Rn/8/V^2;      % RCM
del_R = (1-D)./D * Rn;
core_len = 8;
Srd_RCMC = zeros(Na, Nr);
for i = 1:Na
  for j = 1:Nr
    ideal_p = j + 2*del_R(i,j)/c * fs;          % 需要插值的浮点位置
    quant_p = fix(ideal_p);                         % 需要插值的量化位置
    if (quant_p < core_len/2 || quant_p > Nr - core_len/2)
      Srd_RCMC(i,j) = Srd(i,j);
      continue;
    end
    sig_ps = quant_p + 1 + (-core_len/2:core_len/2-1);      % 信号上的采样点
    sinc_core = sinc(ideal_p - sig_ps);               % 插值核系数
    sinc_core = sinc_core/sum(sinc_core);             % 核系数归一化
    Srd_RCMC(i,j) = sum(sinc_core .* Srd(i, sig_ps)); % 加权求和
  end
end

figure
imagesc(abs(Srd_RCMC))

%Srd_RCMC = Srd;

% *********** 方位压缩 ****************** %
Ha = exp(1j*4*pi/lamda * D * Rn);      % 方位压缩滤波器
Srd2 = Srd_RCMC .* Ha;
Sac = ifft(Srd2, [], 1);                 % 逆变换得二维压缩图像域结果
Sac = Sac(1:na, :);

figure
imagesc(abs(Sac));

figure
contour(abs(Sac))


temp = abs(Sac(fix(na/2), :));
temp = 20*log10(temp/max(temp));

figure
plot(temp)

















end