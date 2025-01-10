function [TransSignal , TxData, Txconst] = SignalGeneration(PhyPar , CPrate , Modulation, WaveShape)
    TxData = randi([0,Modulation.order-1], PhyPar.NumberofCarrier , PhyPar.NumberofSymbol);
    ModData = feval(Modulation.scheme , TxData, Modulation.order, 'Gray', 'Unitaveragepower', true);
    Txconst = ModData;
    CPmat = [zeros(round(CPrate*PhyPar.NumberofCarrier),round((1-CPrate)*PhyPar.NumberofCarrier)), ...
        eye(round(CPrate*PhyPar.NumberofCarrier)); eye(PhyPar.NumberofCarrier)];
    TransSignal = CPmat*sqrt(PhyPar.NumberofCarrier)*WaveShape*ifft(ModData,PhyPar.NumberofCarrier,1);
end

