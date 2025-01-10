clear;close all
% 参数初始化
allpath = genpath(pwd);
addpath(allpath);
[PhyPar,SysPar] = ParameterInitialization(2.8e9, 160, 70, 32, 32, 1, 1, 1, 1, 50, 0.5, 0.5, 2.4e5);
Modulation = struct('scheme', 'qammod' , 'order' , 16); CPrate = 1/16; RCS = 1;
BDPar = struct('waveform',ones((1+CPrate)*PhyPar.NumberofCarrier), ...
    'order',2, 'processdelay', rand(1)*(CPrate)/PhyPar.Subcarrierspacing);
WaveShape = eye(PhyPar.NumberofCarrier);
% 笛卡尔坐标生成
Position = PositionGeneration(PhyPar,SysPar,[15,0,15]);
% 地面镜像
Position.targetground = [1,1,-1].*Position.target;
Position.BDground = [1,1,-1].*Position.target;
Position.userground = [1,1,-1].*Position.user;
Position.scatterground = [1,1,-1].*Position.scatter;

% -------------------------------
%     生成信道
% -------------------------------

pathfirstorder(1,:) = {{'antenna',':'},{'target',1},{'antenna',':'}};
pathfirstorder(2,:) = {{'antenna',':'},{'targetground',1},{'antenna',':'}};

chnfirstord(1) = CHL(Position, PhyPar, pathfirstorder(1,:), RCS);
chnfirstord(2) = CHL(Position, PhyPar, pathfirstorder(2,:), RCS);

% -------------------------------
%     生成信号
% -------------------------------

[TransSignal , TxData, Txconst] = SignalGeneration(PhyPar , CPrate , Modulation,WaveShape);

Delayoffgrid = chnfirstord(1).channel.timeofflight-ceil(chnfirstord(1).channel.timeofflight*PhyPar.NumberofCarrier*PhyPar.Subcarrierspacing)/...
    PhyPar.NumberofCarrier/PhyPar.Subcarrierspacing;
Delayongrid = floor(chnfirstord(1).channel.timeofflight*PhyPar.NumberofCarrier*PhyPar.Subcarrierspacing);

TimeDelaySignal = ReceiveSignalGenerate(Txconst, PhyPar , CPrate , WaveShape, Delayoffgrid(1));
TimeDelaySignal = [zeros(Delayongrid(1),1);reshape(TimeDelaySignal,[numel(TimeDelaySignal),1])];
TransBeamforming = [1;zeros(PhyPar.NumberofAntennahorizon*PhyPar.NumberofAntennavertical-1,1)];
RecSignal = chnfirstord(1).channel.pathloss(PhyPar.NumberofAntennahorizon*...
    PhyPar.NumberofAntennavertical/2)*chnfirstord(1).channel.spatial*TransBeamforming*...
    [zeros(1,Delayongrid(1)),TimeDelaySignal.'];

% -------------------------------
%     空域估计
% -------------------------------

Resolution = [500,500,500];
Num_d =1;
RangeLimit = [0.62*sqrt(PhyPar.Spacinghorizon^3/PhyPar.lambda), ...
    (PhyPar.Spacinghorizon^2*PhyPar.NumberofAntennahorizon^2+...
    PhyPar.Spacingvertical^2*PhyPar.NumberofAntennavertical^2)/PhyPar.lambda];
est = MUSIC2D(PhyPar, Position, RecSignal, RangeLimit,Num_d,Resolution);