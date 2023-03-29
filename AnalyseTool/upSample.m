%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Title: 
% Desc: 
% Author: Jc Chen
% Modified: 2023/03/29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output =  upSample(sig, mul_times)
  [row, col] = size(sig);
  % 一维信号
  if (row == 1 || col == 1)
    sig_len = length(sig);
    S_spec = fft(sig(:));       % 暂时视为列向量
    % 频域补零--时域内插
    S_expand  = [S_spec(1:ceil(sig_len/2)); zeros((mul_times-1)*sig_len, 1); S_spec(ceil(sig_len/2)+1:end)];
    output = ifft(S_expand) * mul_times;         % 调整傅里叶变换的幅度因子
    if (row == 1) 
      output = output.'; 
    end
  
  % 二维信号
  else 
    if (length(mul_times) == 1)
      mul_times = [mul_times, mul_times];
    end
    mul_a = mul_times(1); mul_r = mul_times(2);
    % 直接转换至二维频域
    S_ra = fft(fft(sig, [], 2), [], 1);
    % 构建升采样矩阵
    S_expand = [S_ra(:, 1:ceil(col/2)), zeros(row, (mul_r-1)*col), S_ra(:,ceil(col/2)+1:end)];
    S_expand = [S_expand(1:ceil(row/2),:); zeros((mul_a-1)*row, mul_r*col); S_expand(ceil(row/2)+1,:)];
    output = ifft(ifft(S_expand,[],2), [], 1) * mul_a * mul_r;
  end

end