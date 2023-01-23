test=Clim(stepped,spouseErrH,spouseErrL);
test.M=3;
test=indx_change(test,17);
test=lin_reg(test);
%test.bootstrap(500);
%predicted_test=(Vpn(test.p,test.xL)-test.g2)/log(2.^2);

%X=test.d_mat;
%H=X*inv(X'*test.C_inv*X)*X'*test.C_inv
chist=test.RRbootstrap(500);


