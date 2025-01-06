function [a, b] = factor(N)
    % 检查输入是否为正整数
    if ~isscalar(N) || N <= 0 || floor(N) ~= N
        error('输入必须是一个正整数');
    end

    % 从 sqrt(N) 开始向下搜索因数
    a = floor(sqrt(N));
    while mod(N, a) ~= 0
        a = a - 1;
    end

    % 计算另一个因数
    b = N / a;
end