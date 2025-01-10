function RxSignal = ReceiveSignalGenerate(Txconst, PhyPar , CPrate , WaveShape, Pathdelay)
%RECEIVESIGNALGENERATE 此处显示有关此函数的摘要
%   此处显示详细说明
    CPmat = [zeros(round(CPrate*PhyPar.NumberofCarrier),round((1-CPrate)*PhyPar.NumberofCarrier)), ...
        eye(round(CPrate*PhyPar.NumberofCarrier)); eye(PhyPar.NumberofCarrier)];
    Delaymat = diag(exp(1j*2*pi*(0:PhyPar.NumberofCarrier-1)*(1/PhyPar.NumberofCarrier-Pathdelay*PhyPar.Subcarrierspacing))); 
    RxSignal = CPmat*sqrt(PhyPar.NumberofCarrier)*WaveShape*ifft(Delaymat*Txconst,PhyPar.NumberofCarrier,1);
end

