function position = PositionGeneration(PhyPar, SysPar, BSposition)
    % 参数设置
    NumberofTarget = SysPar.NumberofTarget; % 随机目标点数量
    NumberofScatter = SysPar.NumberofScatter;
    NumberofUser = SysPar.NumberofUser;
    NumberofBD = SysPar.NumberofBD;
    position.BS = BSposition; % UPA阵列的中心坐标
    range = SysPar.Range;
    NumberofAntennahorizon = PhyPar.NumberofAntennahorizon;
    NumberofAntennavertical = PhyPar.NumberofAntennavertical;
    Spacinghorizon = PhyPar.Spacinghorizon;
    Spacingvertical = PhyPar.Spacingvertical;
    
    % UPA阵列参数
    
    position.target = rand(NumberofTarget, 3); 
    position.target = position.target/norm(position.target)*rand(1)*range;
    % -------------------------------
    %     随机生成目标点的坐标
    % -------------------------------

    position.BD = rand(NumberofBD, 3); 
    position.BD = position.BD/norm(position.BD)*rand(1)*range;
    % -------------------------------
    %     随机生成目标点的坐标
    % -------------------------------
    position.scatter = rand(NumberofScatter, 3); 
    position.scatter = position.scatter/norm(position.scatter)*rand(1)*range; 
    % -------------------------------
    %     随机生成用户的坐标
    % -------------------------------
    position.user = rand(NumberofUser, 3); 
    position.user = position.user/norm(position.user)*rand(1)*range; 
    
    % -------------------------------
    %     生成UPA阵列的天线坐标
    % -------------------------------
    % 水平方向和垂直方向的天线索引
    [x_idx, z_idx] = meshgrid(0:NumberofAntennahorizon-1, 0:NumberofAntennavertical-1);
    
    % 计算每个天线单元的相对位置
    x_positions = (x_idx(:) - (NumberofAntennahorizon-1)/2) * Spacinghorizon; % 水平位置
    z_positions = (z_idx(:) - (NumberofAntennavertical-1)/2) * Spacingvertical; % 垂直位置
    
    % UPA阵列坐标（相对于PositionBS）
    position.antenna = [x_positions, zeros(size(x_positions)), z_positions] + position.BS;
end

