function [peak_indices, peak_values] = findPeaks2D(P_dft_angle, N)
    % ȷ�������Ƕ�ά����
    if ndims(P_dft_angle) ~= 2
        error('Input must be a 2D matrix.');
    end
    
    % ��ȡ����Ĵ�С
    [rows, cols] = size(P_dft_angle);
    
    % ����һ���߼��������洢��ֵ
    is_peak = true(rows, cols);
    
    % �������򣨻������ڣ�
    for i = -1:1
        for j = -1:1
            if i == 0 && j == 0
                continue; % �������ĵ�
            end
            % ����������ڵ�ǰ����ƽ��
            shifted_matrix = circshift(P_dft_angle, [i, j]);
            % ��ǰ����Ҫ�������ж�Ӧ���
            is_peak = is_peak & (P_dft_angle > shifted_matrix);
        end
    end
    
    % ��ȡ��ֵ��������Ͷ�Ӧֵ
    [peak_row_indices, peak_col_indices] = find(is_peak);
    peak_values = P_dft_angle(is_peak);
    
    % ���û�з�ֵ��ֱ�ӷ��ؿս��
    if isempty(peak_values)
        peak_indices = [];
        peak_values = [];
        return;
    end
    
    % ����ֵ�Ӵ�С����
    [sorted_values, sort_indices] = sort(peak_values, 'descend');
    sorted_row_indices = peak_row_indices(sort_indices);
    sorted_col_indices = peak_col_indices(sort_indices);
    
    % ѡ��ǰ N ����ֵ������������� N �򷵻����У�
    num_peaks = min(N, numel(sorted_values));
    peak_values = sorted_values(1:num_peaks);
    peak_indices = [sorted_row_indices(1:num_peaks), sorted_col_indices(1:num_peaks)];
end