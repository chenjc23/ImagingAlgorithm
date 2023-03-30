%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Title: 
% Desc: 
% Author: Jc Chen
% Modified: 2023/03/29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function img = RDA(data, Kr, f0, fs, PRF, V, fnc, rng_start, options)
arguments
  data (:,:)
  Kr double
  f0 double
  fs double
  PRF double
  V double
  fnc double
  rng_start double
  options.Tp = nan
  options.theta_bw = nan
  options.type {mustBeMember(options.type, ["point", "scence"])} = "point"
end
Tp = options.Tp;
theta_bw = options.theta_bw;

c = 3e8;
[na, nr] = size(data);
rng_len = 1/fs*c/2 * (nr-1);
lamda = c / f0;
% ************** 二维频域进行一次二次合并距离压缩 ************** %

% 方位向按一个合成孔径时间补零，距离向按一个脉冲长度补零
if (~isnan(theta_bw))
  Ls = theta_bw * (rng_start+rng_len);        % 合成孔径长度
  Ta = Ls / V;                                % 合成孔径时间
  Na = fix(1.2*Ta*PRF) + na;
else
  Na = 2 * na;
end
if (~isnan(Tp))
  tr = rng_start*2/c : 1/fs : (rng_start+rng_len)*2/c + 1.2*Tp;
else
  tr = rng_start*2/c : 1/fs : (rng_start+2*rng_len)*2/c;
end
Nr = length(tr);

Saf = fft(data, Nr, 2);                   % 距离频域
Sdf = fft(Saf, Na, 1);                    % 二维频域

% 构建频率轴                  
fr = fs/Nr * ((0:Nr-1)-fix(Nr/2));            % 距离向频率序列
fr = circshift(fr, -fix(Nr/2));            % 零频移到两端
fa = PRF/Na * ((0:Na-1)-fix(Na/2)).';         
fa = circshift(fa, -fix(Na/2));      
fa = circshift(fa, fix(fnc*Na/PRF)) + fnc;     % 按多普勒中心进行循环移位构建多普勒序列

% 距离压缩
R_ref = rng_start + rng_len/2;                 % 测绘带中心参考距离
D = sqrt(1-lamda^2*fa.^2/4/V^2);               % 徙动因子
Ksrc = 2*V^2*f0^3*D.^3 ./ (c*R_ref*fa.^2);     % 二次调频率
Hr = exp(1j*pi*(1/Kr-1./Ksrc) * fr.^2);        % 一次二次距离压缩滤波器
Srd = ifft(Sdf .* Hr, [], 2);                  % 距离压缩后变换到rd域

% 若为点目标成像则舍去弃置区
if (options.type == "point")
  Srd = Srd(:,1:nr);  
  tr = tr(1:nr);
  Nr = nr;
end

% ************************ RCMC ********************* %
Rn = c/2 * tr;            % 距离序列    
del_R = (1-D)./D * Rn;       % RCM
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

% ****************** 方位压缩 ****************** %
Ha = exp(1j*4*pi/lamda * (D * Rn - repmat(Rn, Na,1)));      % 方位压缩滤波器
Srd2 = Srd_RCMC .* Ha;
img = ifft(Srd2, [], 1);                 % 逆变换得二维压缩图像域结果

% 若为点目标成像，方位向舍弃弃置区得到最终结果
if (options.type == "point")
  img = img(1:na, :);
end

end