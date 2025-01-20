clear;close all
str_vec = @(r) exp(1j*2*pi*r);
L = 2; K_a_half = 16;K_e_half = 16;
K_a = K_a_half*2+1; K_e = K_e_half*2+1;
fc = 2.8e9; c = 3e8; lambda = c/fc;
d_a = 1/2*lambda; d_e = 1/2*lambda;diff_tar_rec = zeros(K_a*K_e,3,L);
azi = zeros(K_a*K_e,L);ele = zeros(K_a*K_e,L);
diff_dis = zeros(K_a*K_e,L);Ant_res = zeros(K_a*K_e,L);

rcs_power = 10; P_noise = -100;% dbm

ant_pos = [kron((-K_a_half:K_a_half)'*d_a,ones(K_e,1)), zeros(K_a*K_e,1), ...
    kron(ones(K_a,1),(-K_e_half:K_e_half)'*d_e)];
tar_pos = [4,10,-5; 1,10,2; 20,-10,2;10,20,1; 14,30,2;];

for i = 1:L
    diff_tar_rec(:,:,i) = tar_pos(i,:)-ant_pos;
    azi(:,i) = atan(diff_tar_rec(:,2,i)./diff_tar_rec(:,1,i));
    ele(:,i) = atan(diff_tar_rec(:,3,i)./sqrt((diff_tar_rec(:,2,i).^2+diff_tar_rec(:,1,i).^2)));
    diff_dis(:,i) = sqrt(sum(diff_tar_rec(:,:,i).^2,2));
    Ant_res(:,i) = str_vec((diff_dis(:,i))/c*fc);
end

N = 2000;
x = randi([0,15],N,L);
s = qammod(x,16,'gray','UnitAveragePower',true);
sgm_rcs = sqrt(rcs_power/2)*(randn(L,1)+1j*randn(L,1));
rec_sig = Ant_res*diag(sgm_rcs)*s.'+sqrt(10^(P_noise/10)/2)*(randn(K_a*K_e,N)+1j*randn(K_a*K_e,N));
expectationFourthOrder = (rec_sig.*conj(rec_sig(end:-1:1,:)))*(rec_sig.*conj(rec_sig(end:-1:1,:)))'/N;
expectationSecondOrder1 = (rec_sig*rec_sig').*(conj(rec_sig(end:-1:1,:))*conj(rec_sig(end:-1:1,:))')/N^2;
templateCount = sum(rec_sig.*conj(rec_sig(end:-1:1,:)),2)/N;
expectationSecondOrder2 = templateCount*templateCount';
Cum_mat = expectationFourthOrder - expectationSecondOrder2 - expectationSecondOrder1;
[Ue,ev] = eig(Cum_mat); 
Num_d = L;
U_sig = Ue(:,1:Num_d);
CHOOSINGMAT1FORWORD = kron([eye(K_a-1),zeros(K_a-1,1)],eye(K_e));
CHOOSINGMAT1BACKWORD= kron([zeros(K_a-1,1),eye(K_a-1)],eye(K_e));
CHOOSINGMAT2FORWORD = kron(eye(K_a),[eye(K_e-1),zeros(K_e-1,1)]);
CHOOSINGMAT2BACKWORD = kron(eye(K_a),[zeros(K_e-1,1),eye(K_e-1)]);
A1Forward = CHOOSINGMAT1FORWORD*U_sig;
A1Backward = CHOOSINGMAT1BACKWORD*U_sig;
A2Forward = CHOOSINGMAT2FORWORD*U_sig;
A2Backward = CHOOSINGMAT2BACKWORD*U_sig;
[~,L1] =  eig((A1Forward'*A1Forward)\A1Forward'*A1Backward);
[~,L2] =  eig((A2Forward'*A2Forward)\A2Forward'*A2Backward);
thetaModify = angle(diag(conj(L1)))./sqrt(4*pi^2-angle(diag(conj(L2))).^2);
thetaEst= acos(thetaModify)/pi*180;
eleEst= acos(angle(diag(conj(L2)))/pi/2)/pi*180;
estimateResult = [thetaEst,eleEst];
estimateResult  = sortrows(estimateResult , -2);
realValue = [azi(K_a_half*K_e+K_e_half+1,:)'/pi*180,90-ele(K_a_half*K_e+K_e_half+1,:)'/pi*180];
realValue  = sortrows(realValue , -2);
disp('Theta Estimation, Theta Real Number, Element Estimation, Element Real Value')
disp([estimateResult(:,1),realValue(:,1),estimateResult(:,2),realValue(:,2)])

[Ue,ev] = eig(Cum_mat); 
Num_d = L;
U_sig = Ue(:,1:Num_d);
CHOOSINGMAT1FORWORD = kron([eye(K_a-1),zeros(K_a-1,1)],eye(K_e));
CHOOSINGMAT1BACKWORD= kron([zeros(K_a-1,1),eye(K_a-1)],eye(K_e));
CHOOSINGMAT2FORWORD = kron(eye(K_a),[eye(K_e-1),zeros(K_e-1,1)]);
CHOOSINGMAT2BACKWORD = kron(eye(K_a),[zeros(K_e-1,1),eye(K_e-1)]);
A1Forward = CHOOSINGMAT1FORWORD*U_sig;
A1Backward = CHOOSINGMAT1BACKWORD*U_sig;
A2Forward = CHOOSINGMAT2FORWORD*U_sig;
A2Backward = CHOOSINGMAT2BACKWORD*U_sig;
[~,L1] =  eig((A1Forward'*A1Forward)\A1Forward'*A1Backward);
[~,L2] =  eig((A2Forward'*A2Forward)\A2Forward'*A2Backward);
thetaModify = angle(diag(conj(L1)))./sqrt(4*pi^2-angle(diag(conj(L2))).^2);
thetaEst= acos(thetaModify)/pi*180;
eleEst= acos(angle(diag(conj(L2)))/pi/2)/pi*180;
estimateResult = [thetaEst,eleEst];
estimateResult  = sortrows(estimateResult , -2);
realValue = [azi((K_a_half)*K_e+K_e_half+1,:)'/pi*180,90-ele((K_a_half)*K_e+K_e_half+1,:)'/pi*180];
realValue  = sortrows(realValue , -2);
disp('Theta Estimation, Theta Real Number, Element Estimation, Element Real Value')
disp([estimateResult(:,1),realValue(:,1),estimateResult(:,2),realValue(:,2)])