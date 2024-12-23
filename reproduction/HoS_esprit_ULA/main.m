clear;close all
str_vec = @(r) exp(1j*2*pi*r);
L = 5;
K_a_half = 64;
K_a = K_a_half*2+1; K_e = 1;
fc = 2.8e10;
c = 3e8;
lambda = c/fc;
d_a = 1/2*lambda; d_e = 1/2*lambda;
diff_tar_rec = zeros(K_a*K_e,3,L);
azi = zeros(K_a*K_e,L);
ele = zeros(K_a*K_e,L);
diff_dis = zeros(K_a*K_e,L);
Ant_res = zeros(K_a*K_e,L);

rcs_power = 10;
P_noise = 0;

ant_pos = [zeros(K_a*K_e,1),kron((0:K_a-1)'*d_a,ones(K_e,1)), kron(ones(K_a,1),(0:K_e-1)'*d_e)];
tar_pos = [30,15,3;...
    10,-2,2;
    20,-15,2;
    10,20,1;
    15,0,0];
for i = 1:L
    diff_tar_rec(:,:,i) = tar_pos(i,:)-ant_pos;
    azi(:,i) = atan(diff_tar_rec(:,2,i)./diff_tar_rec(:,1,i));
    ele(:,i) = atan(diff_tar_rec(:,3,i)./sqrt((diff_tar_rec(:,2,i).^2+diff_tar_rec(:,1,i).^2)));
    diff_dis(:,i) = sqrt(sum(diff_tar_rec(:,:,i).^2,2));
    Ant_res(:,i) = str_vec((diff_dis(:,i)-min(diff_dis(K_a_half+1,i)))/c*fc);
end

N = 512;
x = randi([0,3],N,L);
s = qammod(x,4,'gray','UnitAveragePower',true);
sgm_rcs = sqrt(rcs_power/2)*(randn(L,1)+1j*randn(L,1));
rec_sig = Ant_res*diag(sgm_rcs)*s.'+sqrt(10^(P_noise/10)/2)*(randn(K_a*K_e,N)+1j*randn(K_a*K_e,N));

[Ue,ev] = eig(rec_sig * rec_sig');
Num_d = L;
U_noise = Ue(:,1:K_a*K_e-Num_d);
U_sig = Ue(:,(K_a*K_e-Num_d+1):K_a*K_e);

%% estimation
% grid
far_limit = d_a^2*K_a^2/lambda;
sense_limit = 0.62*sqrt(d_a^3/lambda);
r = (0:1/500:1)*(far_limit-sense_limit)+sense_limit;
theta = (-0.5:1/400:0.5)*2;
test_range = -(-K_a_half:K_a_half)'.*permute(theta'*ones(size(r)),[3,1,2])/2 +...
     d_a^2/2/lambda*((-K_a_half:K_a_half).^2)'.*permute((1-theta.^2)'./r,[3,1,2]);
t_bas_2nd = str_vec(test_range);
P_mu = zeros(401,501);
for i = 1:401
    for j = 1:501
        P_mu(i,j) = (t_bas_2nd(:,i,j)' * (U_noise*U_noise') * t_bas_2nd(:,i,j));
    end
end
P_pod = 1./P_mu;
[a,pc] = polarPcolor(r,asin(theta)/pi*180,100*abs(P_pod)/max(max(abs(P_pod))),'colBar',0);
scatter(sin(azi(65,:)).*diff_dis(65,:)/far_limit,cos(azi(65,:)).*diff_dis(65,:)/far_limit, 'r*');
r_est = zeros(Num_d,1);
theta_est = zeros(Num_d,1);
for count = 1:Num_d
    P_pod = 1./P_mu;
    [~,idx]=max(P_pod(:));
    r_est(count) = floor(idx/401)/500*(far_limit-sense_limit)+sense_limit;
    stheta = (mod(idx-1, 401)-200)/200;
    theta_est(count) = asin((mod(idx-1, 401)-200)/200);
    v_recon = str_vec(-(-K_a_half:K_a_half)'*stheta/2 +...
    d_a^2/2/lambda*((-K_a_half:K_a_half).^2)'*(1-stheta^2)'/r_est(count));
    P_mu = P_mu+abs(permute(pagemtimes(permute(conj(t_bas_2nd),[2,1,3]),v_recon),[1,3,2])).^2;

end
scatter(r_est.*sin(theta_est)/far_limit,r_est.*cos(theta_est)/far_limit, 'go');
% %% DFT
% P_dft = zeros(401,501);
% for i = 1:401
%     for j = 1:501
%         P_dft(i,j) = sqrt(sum(abs(t_bas_2nd(:,i,j)'*rec_sig).^2,2));
%     end
% end
% figure()
% [a,pc] = polarPcolor(r,asin(theta)/pi*180,abs(P_dft)/max(max(abs(P_dft))));
% title('DFT')


