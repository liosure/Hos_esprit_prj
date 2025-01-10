function estidx = MUSIC2D(PhyPar, Position, signal, rangelimit,Num_d,resolution)
% estimation
% grid
K_e = PhyPar.NumberofAntennavertical;
K_a = PhyPar.NumberofAntennahorizon;
lambda = PhyPar.lambda;
sense_limit = rangelimit(1);
far_limit = rangelimit(2);
[Ue,~] = eig(signal * signal');
U_noise = Ue(:,1:K_a*K_e-Num_d);
r_res = resolution(1);
theta_res = resolution(2);
phi_res = resolution(3);
r = (0:1/r_res:1-1/r_res)'*(far_limit-sense_limit)+sense_limit;
theta = (-0.5:1/theta_res:0.5-1/theta_res)'*2;
phi = (-0.5:1/phi_res:0.5-1/phi_res)'*2;
Ant_pos = Position.antenna-Position.BS;
TestRangefun = @(theta, phi, r)  -phi * Ant_pos(:,3)'/lambda - sqrt(1-phi.^2) .* theta * Ant_pos(:,1)'/lambda +  ...
     (1-phi.^2)./r * Ant_pos(:,3)'.^2/2/lambda + (1-phi.^2).*theta.^2./r * Ant_pos(:,1)'.^2/2/lambda - ...
     -phi .* sqrt(1-phi.^2) .* theta./r*(Ant_pos(:,1).*Ant_pos(:,3))'/lambda;
estidx = [round(theta_res/2), round(phi_res/2), round(r_res/2)];
estidxorigin = [1,1,1];
while sum(estidx-estidxorigin) 
    estidxorigin = estidx;
    test_range = TestRangefun(theta,phi(estidx(2)),r(estidx(3)));
    testbasic = PhyPar.StrVecfun(test_range);
    P_mu = 1./abs(sum(abs( testbasic * U_noise).^2,2));
    [~,estidx(1)] = max(P_mu);
    test_range = TestRangefun(theta(estidx(1)),phi,r(estidx(3)));
    testbasic = PhyPar.StrVecfun(test_range);
    P_mu = 1./abs(sum(abs( testbasic * U_noise).^2,2));
    [~,estidx(2)] = max(P_mu);
    test_range = TestRangefun(theta(estidx(1)),phi(estidx(2)),r);
    testbasic = PhyPar.StrVecfun(test_range);
    P_mu = 1./abs(sum(abs( testbasic * U_noise).^2,2));
    [~,estidx(3)] = max(P_mu);
end

end

