function sra = getSimulateEcho(azm_len, rng_len, rng_start, V, PRF, theta_rc, La, f0, Tp, B, fs, tg_pos)

c = 3e8;
lamda = c/fc;
theta_bw = 0.886*lamda/La;
Kr = B/Tp;


tr = rng_start*2/c : 1/fs : (rng_start+rng_len)*2/c;          % 距离向时间序列
Nr = length(tr);
Na = fix(azm_len/V * PRF / 2) * 2 + 1;                        % 方位向采样点数（保证能构成关于中点对称）
Na_vec = linspace(-(Na-1)/2, (Na-1)/2, Na);

sra = zeros(Na, Nr);              % 存储回波信号
for n = 1:Na
  % 雷达坐标，默认雷达距离向坐标为0
  radar_x = 0; 
  radar_y = V/PRF * Na_vec(n);
  % 对每条距离线遍历所有的目标
  for i = 1:size(tg_pos, 1)
    tg_x = tg_pos(i, 1); tg_y = tg_pos(i, 2);                    % 获取目标坐标                           
    R = sqrt((tg_x-radar_x)^2 + (tg_y-radar_y)^2);               % 雷达与目标距离
    theta = atan((tg_y-radar_y)/(tg_x-radar_x)) - theta_rc;     % 目标在斜距平面与视线的夹角
    wa = sinc(0.886*theta/theta_bw)^2;                            % 方位包络（某一方位时刻为常数）
    wr = (tr-2*R/c) <= Tp;                                       % 距离包络（构造窗序列）
    sra(n,:) = sra(n,:) + wa * wr .* exp(-1j*4*pi*f0*R/c) .* exp(1j*pi*Kr*(tr-2*R/c).^2);        % 当前目标的回波
  end
end

