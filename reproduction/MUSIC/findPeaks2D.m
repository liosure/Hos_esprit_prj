function [peak_indices, peak_values] = findPeaks2D(P_dft_angle, N)
    % 确保输入是二维矩阵
    if ndims(P_dft_angle) ~= 2
        error('Input must be a 2D matrix.');
    end
    
    % 获取矩阵的大小
    [rows, cols] = size(P_dft_angle);
    
    % 创建一个逻辑矩阵来存储峰值
    is_peak = true(rows, cols);
    
    % 遍历邻域（滑动窗口）
    for i = -1:1
        for j = -1:1
            if i == 0 && j == 0
                continue; % 跳过中心点
            end
            % 将矩阵相对于当前方向平移
            shifted_matrix = circshift(P_dft_angle, [i, j]);
            % 当前点需要比邻域中对应点大
            is_peak = is_peak & (P_dft_angle > shifted_matrix);
        end
    end
    
    % 提取峰值点的索引和对应值
    [peak_row_indices, peak_col_indices] = find(is_peak);
    peak_values = P_dft_angle(is_peak);
    
    % 如果没有峰值，直接返回空结果
    if isempty(peak_values)
        peak_indices = [];
        peak_values = [];
        return;
    end
    
    % 按峰值从大到小排序
    [sorted_values, sort_indices] = sort(peak_values, 'descend');
    sorted_row_indices = peak_row_indices(sort_indices);
    sorted_col_indices = peak_col_indices(sort_indices);
    
    % 选择前 N 个峰值（如果数量不足 N 则返回所有）
    num_peaks = min(N, numel(sorted_values));
    peak_values = sorted_values(1:num_peaks);
    peak_indices = [sorted_row_indices(1:num_peaks), sorted_col_indices(1:num_peaks)];
end