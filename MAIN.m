clear;close all
% 参数初始化
allpath = genpath(pwd);
addpath(allpath);
[physicalParameter,systemParameter] = ParameterInitialization(2.8e9, 160, ...
    70, 32, 32, 1, 1, 1, 1, 50, 0.5, 0.5, 2.4e5);
% 参数说明：
modulationParameter = struct('scheme', 'qammod' , 'order' , 16); CPrate = 1/16; RCS = 1;
BDParameter = struct('waveform',ones((1+CPrate)*physicalParameter.NumberofCarrier), ...
    'order',2, 'processdelay', rand(1)*(CPrate)/physicalParameter.Subcarrierspacing);
waveShape = eye(physicalParameter.NumberofCarrier);
% 笛卡尔坐标生成 Freqc,NumberofCarrier,NumberofSymbol, ...
% NumberofAntennahorizon, NumberofAntennavertical, NumberofTarget, NumberofScatter, NumberofBD, ...
% NumberofUser, Range, HorizonSpacingRate, VerticalSpacingRate , Subcarrierspacing
POSITION = PositionGeneration(physicalParameter,systemParameter,[15,0,15]);
% 地面镜像
POSITION.targetground = [1,1,-1].*POSITION.target;
POSITION.BDground = [1,1,-1].*POSITION.target;
POSITION.userground = [1,1,-1].*POSITION.user;
POSITION.scatterground = [1,1,-1].*POSITION.scatter;
% -------------------------------
%     生成信道
% -------------------------------
pathFirstOrder(1,:) = {{'antenna',':'},{'target',1},{'antenna',':'}};
pathFirstOrder(2,:) = {{'antenna',':'},{'targetground',1},{'antenna',':'}};
channelFirstOrder(1) = CHL(POSITION, physicalParameter, pathFirstOrder(1,:), RCS);
channelFirstOrder(2) = CHL(POSITION, physicalParameter, pathFirstOrder(2,:), RCS);
% -------------------------------
%     生成发射信号
% -------------------------------
[txSignal , txData, txConst] = SignalGeneration(physicalParameter, ...
    CPrate, modulationParameter, waveShape);
% -------------------------------
%     生成反射接收信号
% -------------------------------
sampleNumber = (1+CPrate)*physicalParameter.NumberofCarrier * (physicalParameter.NumberofSymbol+1);
delayOffGrid = channelFirstOrder(1).channel.timeofflight-...
    ceil(channelFirstOrder(1).channel.timeofflight*...
    physicalParameter.NumberofCarrier*physicalParameter.Subcarrierspacing)/...
    physicalParameter.NumberofCarrier/physicalParameter.Subcarrierspacing;
delayOnGrid = floor(channelFirstOrder(1).channel.timeofflight*...
    physicalParameter.NumberofCarrier*physicalParameter.Subcarrierspacing);

timeDelaySignal = ReceiveSignalGenerate(txConst, physicalParameter ,...
    CPrate , waveShape, delayOffGrid(1));
timeDelaySignal = [zeros(delayOnGrid(1),1);reshape(timeDelaySignal,[numel(timeDelaySignal),1])];
txBeamforming = [1;zeros(physicalParameter.NumberofAntennahorizon*...
    physicalParameter.NumberofAntennavertical-1,1)];
rxSignal = channelFirstOrder(1).channel.pathloss(physicalParameter.NumberofAntennahorizon*...
    physicalParameter.NumberofAntennavertical/2)*channelFirstOrder(1).channel.spatial*...
    txBeamforming*[zeros(1,delayOnGrid(1)),timeDelaySignal.',...
    zeros(1,sampleNumber-delayOnGrid(1)-numel(timeDelaySignal))];

% -------------------------------
%     空域估计
% -------------------------------

% RESLUTION = [500,500,500];
% detectTargetNumber =1;
% rangeLimit = [0.62*sqrt(physicalParameter.Spacinghorizon^3/physicalParameter.lambda), ...
%     (physicalParameter.Spacinghorizon^2*physicalParameter.NumberofAntennahorizon^2+...
%     physicalParameter.Spacingvertical^2*physicalParameter.NumberofAntennavertical^2)/...
%     physicalParameter.lambda];
% estimate = MUSIC2D(physicalParameter, POSITION, rxSignal, ...
%     rangeLimit,detectTargetNumber,RESLUTION);