%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Title: 
% Desc: 
% Author: Jc Chen
% Modified: 2023/03/29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function inspectView(sig, upsample_times, win_len)
arguments
  sig (:,:)
  upsample_times (:,:)
  win_len (:,:) = [64, 64]
end

  [row, col] = size(sig);
  sig_abs = abs(sig);
  [~, loc_max] = max(sig_abs(:));
  loc_r = ceil(loc_max / row);
  loc_c = mod(loc_max, row);
  
  % 对信号进行窗截取
  row_start = loc_r - win_len(1)/2 + 1;
  row_end   = loc_r + win_len(1)/2;
  col_start = loc_c - win_len(2)/2 + 1;
  col_end   = loc_c + win_len(2)/2;
  if (row_start < 1) 
    row_start = 1; 
  end
  if (row_end > row) 
    row_end = row; 
  end
  if (col_start < 1)
    col_start = 1;
  end
  if (col_end > col)
    col_end = col;
  end
  sig_win = sig(row_start:row_end, col_start:col_end);
  
  % 信号升采样
  sig_expand = upSample(sig_win, upsample_times);
  
  % 归一化作图
  sig_expand = abs(sig_expand);
  
  if (row == 1 || col == 1)
    sig_expand = 20*log10(sig_expand/max(sig_expand(:)));
    figure;
    plot(sig_expand, 'Color', [1 0.5 0], LineWidth=1); grid on;
    axis tight;
  else
    figure;
    contour(sig_expand, 20);
    title("等高线图像");
  
    figure;
    image(255*sig_expand/max(sig_expand(:)));
    colormap("gray");
    title("灰度图");
  end
end