classdef BDCHL < CHL
    %UNTITLED 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        BDsignal
        BDsymbol
    end
    
    methods
        function obj = BDCHL(position, PhyPar, path, RCS, BDPar)
            %UNTITLED 构造此类的实例
            %   此处显示详细说明
            obj@CHL(position, PhyPar, path, RCS);
            obj.ProcessDelay = ProcDelay;
            obj = obj.BDSignalGeneration(obj, PhyPar, BDPar);
        end
        
        function obj = BDSignalGeneration(obj, PhyPar, BDPar)
            %METHOD1 此处显示有关此方法的摘要
            %   此处显示详细说明
            obj.BDsymbol = randi([0,BDPar.order-1],PhyPar.NumberofSymbol,1);
            signal = BDPar.waveform(:,obj.BDsymbol+1);
            obj.BDsignal = [zeros(floor((BDPar.processdelay+obj.firstToF)/(PhyPar.Subcarrierspacing*PhyPar.NumberofCarrier)),1);...
                reshape(signal, [numel(signal),1])];
        end
    end
end

