classdef CHL
    %CHL 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        dep
        midpath
        arr
        pos
        channel
        firstToF
    end
    
    methods
        function obj = CHL(position, PhyPar, path, RCS)
            %CHL 构造此类的实例
            %   此处显示详细说明
            obj.pos = position;
            Numberofpath = numel(path)-1;
            % ChannelStruct = struct('azi', [], 'ele', [], 'timeofflight', [],'pathloss',[],'antennaresponse',[]);
            % obj.channel = repmat(ChannelStruct, Numberofpath, 1);
            tempchannel = obj.chl_gen(path{1},path{2},PhyPar);
            obj.channel.azi(1) = tempchannel.azi(floor(end/2));
            obj.channel.ele(1) = tempchannel.ele(floor(end/2));
            obj.channel.timeofflight = tempchannel.timeofflight;
            obj.firstToF = tempchannel.timeofflight;
            obj.channel.pathloss = tempchannel.pathloss;
            obj.channel.spatial = tempchannel.antennaresponse.';
            for i = 2:Numberofpath
                tempchannel = obj.chl_gen(path{i},path{i+1},PhyPar);
                obj.channel.azi(i) = tempchannel.azi(floor(end/2));
                obj.channel.ele(i) = tempchannel.ele(floor(end/2));
                obj.channel.timeofflight = obj.channel.timeofflight + tempchannel.timeofflight;
                obj.channel.pathloss = obj.channel.pathloss.*tempchannel.pathloss.*RCS(i-1);
                obj.channel.spatial = tempchannel.antennaresponse*obj.channel.spatial;
            end
        end
        
        function Channel = chl_gen(obj,poscell1,poscell2,PhyPar)
            % cell：阵列响应矩阵，水平角，俯仰角，时延
            pos1 = cell2mat(poscell1(1));
            index1 = cell2mat(poscell1(2));
            pos2 = cell2mat(poscell2(1));
            index2 = cell2mat(poscell2(2));
            Posdiff = obj.pos.(pos2)(index2,:)-obj.pos.(pos1)(index1,:);
            Channel.azi = atan(Posdiff(:,2)./Posdiff(:,1));
            Channel.ele = atan(Posdiff(:,3)./sqrt((Posdiff(:,2).^2+Posdiff(:,1).^2)));
            Channel.timeofflight = sqrt(sum(Posdiff.^2,2))/PhyPar.c;
            Channel.pathloss = sqrt(PhyPar.lambda^2/(4*pi)^2./sum(Posdiff.^2,2));
            if strcmp(pos2, 'antenna')
                diff_dis_center_tar = sqrt(sum((obj.pos.BS-obj.pos.(pos1)(index1,:)).^2,2))/PhyPar.c;
                Channel.antennaresponse = PhyPar.StrVecfun((Channel.timeofflight-diff_dis_center_tar)*PhyPar.Freqc);
            elseif strcmp(pos1, 'antenna')
                diff_dis_center_tar = sqrt(sum((obj.pos.(pos2)(index2,:) - obj.pos.BS).^2,2))/PhyPar.c;
                Channel.antennaresponse = PhyPar.StrVecfun((Channel.timeofflight-diff_dis_center_tar)*PhyPar.Freqc);
            else
                Channel.antennaresponse = 1;
            end
        end
    end
end

