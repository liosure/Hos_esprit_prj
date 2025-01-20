clear;close all
str_vec = @(r) exp(1j*2*pi*r);
L = 1; K_a_half = 16;K_e_half = 16;
K_a = K_a_half*2+1; K_e = K_e_half*2+1;
fc = 2.8e9; c = 3e8; lambda = c/fc;
d_a = 1/2*lambda; d_e = 1/2*lambda;diff_tar_rec = zeros(K_a*K_e,3,L);
azi = zeros(K_a*K_e,L);ele = zeros(K_a*K_e,L);
diff_dis = zeros(K_a*K_e,L);Ant_res = zeros(K_a*K_e,L);

rcs_power = 1; P_noise = -200;% dbm

ant_pos = [zeros(K_a*K_e,1),kron((-K_a_half:K_a_half)'*d_a,ones(K_e,1)), kron(ones(K_a,1),(0:K_e-1)'*d_e)];
tar_pos = [10,4,1; 10,-2,2; 20,-10,2;10,20,1; 14,30,2;];

for i = 1:L
    diff_tar_rec(:,:,i) = tar_pos(i,:)-ant_pos;
    azi(:,i) = atan(diff_tar_rec(:,2,i)./diff_tar_rec(:,1,i));
    ele(:,i) = atan(diff_tar_rec(:,3,i)./sqrt((diff_tar_rec(:,2,i).^2+diff_tar_rec(:,1,i).^2)));
    diff_dis(:,i) = sqrt(sum(diff_tar_rec(:,:,i).^2,2));
    Ant_res(:,i) = str_vec((diff_dis(:,i)-min(diff_dis(K_a_half+1,i)))/c*fc);
end

N = 2048;
x = randi([0,15],N,L);
s = qammod(x,16,'gray','UnitAveragePower',true);
sgm_rcs = sqrt(rcs_power/2)*(randn(L,1)+1j*randn(L,1));
rec_sig = Ant_res*diag(sgm_rcs)*s.'+sqrt(10^(P_noise/10)/2)*(randn(K_a*K_e,N)+1j*randn(K_a*K_e,N));

idxm = ((K_a_half+1):(2*K_a_half))-1;
idxm1 = ((K_a_half+2):(2*K_a_half+1))-1;
idxmm1 = ((K_a_half):(2*K_a_half-1))-1;
idxn = ((K_a_half+1):(2*K_a_half))-1;
idxn1 = ((K_a_half+2):(2*K_a_half+1))-1;
idxmn = ((K_a_half+1):-1:2)-1;
idx1mn = ((K_a_half+2):-1:3)-1;


FourthOrderExpectation1 = (conj(rec_sig(idxm,:)).*rec_sig(idxm1,:))*...
    (conj(rec_sig(idxn,:)).*rec_sig(idxn1,:))'/N;
SecondOrderExpectation1 = sum(conj(rec_sig(idxm,:)).*rec_sig(idxm1,:),2)*...
    sum(conj(rec_sig(idxn1,:)).*rec_sig(idxn,:),2).'/N^2;
SecondOrderExpectation2 = (conj(rec_sig(idxm,:))*rec_sig(idxn,:).').*...
    (rec_sig(idxm1,:)*rec_sig(idxn1,:)')/N^2;
SecondOrderExpectation3 = (conj(rec_sig(idxm,:))*rec_sig(idxn1,:)').*...
    (rec_sig(idxm1,:)*rec_sig(idxn,:).')/N^2;

Cum1 = FourthOrderExpectation1-SecondOrderExpectation1-...
    SecondOrderExpectation2-SecondOrderExpectation3;

FourthOrderExpectation2 = (conj(rec_sig(idxmm1,:)).*rec_sig(idxm,:))*...
    (conj(rec_sig(idx1mn,:)).*rec_sig(idxmn,:))'/N;
SecondOrderExpectation1 = sum(conj(rec_sig(idxmm1,:)).*rec_sig(idxm,:),2)*...
    sum(conj(rec_sig(idxmn,:)).*rec_sig(idx1mn,:),2).'/N^2;
SecondOrderExpectation2 = (conj(rec_sig(idxmm1,:))*rec_sig(idx1mn,:).').*...
    (rec_sig(idxm,:)*rec_sig(idxmn,:)')/N^2;
SecondOrderExpectation3 = (conj(rec_sig(idxmm1,:))*rec_sig(idxmn,:)').*...
    (rec_sig(idxm,:)*rec_sig(idx1mn,:).')/N^2;

Cum2 = FourthOrderExpectation2-SecondOrderExpectation1-...
    SecondOrderExpectation2-SecondOrderExpectation3;

FourthOrderExpectation3 = (conj(rec_sig(idxm,:)).*rec_sig(idxm1,:))*...
    (conj(rec_sig(idx1mn,:)).*rec_sig(idxmn,:))'/N;
SecondOrderExpectation1 = sum(conj(rec_sig(idxm,:)).*rec_sig(idxm1,:),2)*...
    sum(conj(rec_sig(idxmn,:)).*rec_sig(idx1mn,:),2).'/N^2;
SecondOrderExpectation2 = (conj(rec_sig(idxm,:))*rec_sig(idx1mn,:).').*...
    (rec_sig(idxm1,:)*rec_sig(idxmn,:)')/N^2;
SecondOrderExpectation3 = (conj(rec_sig(idxm,:))*rec_sig(idxmn,:)').*...
    (rec_sig(idxm1,:)*rec_sig(idx1mn,:).')/N^2;

Cum3 = FourthOrderExpectation3-SecondOrderExpectation1-...
    SecondOrderExpectation2-SecondOrderExpectation3;

FourthOrderExpectation4 = (conj(rec_sig(idxmm1,:)).*rec_sig(idxm,:))*...
    (conj(rec_sig(idxn,:)).*rec_sig(idxn1,:))'/N;
SecondOrderExpectation1 = sum(conj(rec_sig(idxmm1,:)).*rec_sig(idxm,:),2)*...
    sum(conj(rec_sig(idxn1,:)).*rec_sig(idxn,:),2).'/N^2;
SecondOrderExpectation2 = (conj(rec_sig(idxmm1,:))*rec_sig(idxn,:).').*...
    (rec_sig(idxm,:)*rec_sig(idxn1,:)')/N^2;
SecondOrderExpectation3 = (conj(rec_sig(idxmm1,:))*rec_sig(idxn1,:)').*...
    (rec_sig(idxm,:)*rec_sig(idxn,:).')/N^2;

Cum4 = FourthOrderExpectation4-SecondOrderExpectation1-...
    SecondOrderExpectation2-SecondOrderExpectation3;


Cum_mat = [Cum1,Cum4,Cum2;...
    Cum4',Cum1,Cum3;...
    Cum2',Cum3',Cum1];

Exp_mat = [FourthOrderExpectation1,FourthOrderExpectation4,FourthOrderExpectation2;...
    FourthOrderExpectation4',FourthOrderExpectation1,FourthOrderExpectation3;...
    FourthOrderExpectation2',FourthOrderExpectation3',FourthOrderExpectation1];

[Ue,ev] = eig(Cum_mat);
Num_d = L;
U_noise = Ue(:,1:K_a*K_e-Num_d);
U_sig = Ue(:,1:Num_d);

A1 = U_sig(1:K_a_half,:);
A2 = U_sig(K_a_half+1:2*K_a_half,:);
A3 = U_sig(2*K_a_half+1:3*K_a_half,:);

[~,L1] = eig((A1'*A1)\A1'*A2);
[~,L2] = eig((A1'*A1)\A1'*A3);

Thetaest= asin(angle(diag(L2))/pi/2)/pi*180;
disp([Thetaest, azi(K_a_half+1,1)/pi*180])
Rangeest = 1./(angle(diag(L1))/2/pi./cos(Thetaest/180*pi).^2/d_a^2*lambda);



