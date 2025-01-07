clear;close all
SteeringVectorfun = @(r) exp(1j*2*pi*r);
Freqc = 2.8e10; c = 3e8; lambda = c/Freqc;
NumberofCarrier = 16;
NumberofSymbol = 16;


% 参数设置
NumberofTarget = 2; % 随机目标点数量
NumberofScatter = 3;
NumberofUser = 3;
PositionBS = [0, 0, 0]; % UPA阵列的中心坐标

% UPA阵列参数
NumberofAntennahorizon = 4; % 水平方向天线数量
Spacinghorizon = 0.5*lambda; % 水平方向天线间距（单位：波长）
NumberofAntennavertical = 3; % 垂直方向天线数量
Spacingvertical = 0.5*lambda; % 垂直方向天线间距（单位：波长）

% -------------------------------
%     随机生成目标点的坐标
% -------------------------------
TargetPositions = rand(NumberofTarget, 3) .* [50,50,50]; % 目标点在(0,0,0)到(100,100,100)范围内随机生成

% -------------------------------
%     随机生成目标点的坐标
% -------------------------------
ScatterPositions = rand(NumberofScatter, 3) .* [50,50,50]; % 目标点在(0,0,0)到(100,100,100)范围内随机生成

% -------------------------------
%     随机生成用户的坐标
% -------------------------------
UserPositions = rand(NumberofUser, 3) .* [50,50,50]; % 目标点在(0,0,0)到(100,100,100)范围内随机生成

% -------------------------------
%     生成UPA阵列的天线坐标
% -------------------------------
% 水平方向和垂直方向的天线索引
[x_idx, z_idx] = meshgrid(0:NumberofAntennahorizon-1, 0:NumberofAntennavertical-1);

% 计算每个天线单元的相对位置
x_positions = (x_idx(:) - (NumberofAntennahorizon-1)/2) * Spacinghorizon; % 水平位置
z_positions = (z_idx(:) - (NumberofAntennavertical-1)/2) * Spacingvertical; % 垂直位置

% UPA阵列坐标（相对于PositionBS）
Antenna_Positions = [x_positions, zeros(size(x_positions)), z_positions] + PositionBS;

% -------------------------------
%     生成信道
% -------------------------------

Channel = CHL(pos_tar_input, pos_ant_input, pos_ref_input,K_a_half,c,fc,Num_ref);









