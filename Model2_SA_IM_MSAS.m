function [rho50, sigma_z] = Model2_SA_IM_MSAS(X1,X2,X3,X4,X5)

% Created by Mao-Xin Wang (dr.maoxin.wang@gmail.com or maoxinwang@ust.hk)
% January 2025
%
% Model 2: Correlations between SAs and other IMs (SA-IM pairs)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%
%   X1   =  scalar or matrix of T (for SA)
%   X2   =  scalar or matrix of ksi (damping ratio, in percentage)
%   X3   =  scalar or matrix of I_(AS,SA) (set to 1 when SA corresponds to AS; 
%                                          set to 0 when SA corresponds to MS)
%   X4   =  scalar or matrix of I_IM (set to 1, 2, ..., or 9 when IM denotes
%             ASI,SI, DSI, PGA, PGV, AI, CAV, Ds5-75, or Ds5-95, respectively)
%   X5   =  scalar or matrix of I_HV (set to 1, 2, 3, or 4 for 
%           MS_H-AS_H, MS_H-AS_V, MS_V-AS_H,or MS_V-AS_V, respectively)
%          (Note: the above inputs must be in the same dimension)
%
% OUTPUT
%
%   rho50     =  scalar or matrix of median correlation
%   sigma_z   =  scalar or matrix of standard deviation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load coeff_Model2
X_min = [log(0.01) 	log(0.5) 	0 	1 	1];
X_max = [log(5) 	log(30) 	1 	9 	4];
[n_row,n_col] = size(X1);
n_data = n_row*n_col;
X1_temp = reshape(log(X1),n_data,1);
X2_temp = reshape(log(X2),n_data,1);
X3_temp = reshape(X3,n_data,1);
X4_temp = reshape(X4,n_data,1);
X5_temp = reshape(X5,n_data,1);
X_norm = ([X1_temp,X2_temp,X3_temp,X4_temp,X5_temp]-repmat(X_min,[n_data,1]))./(repmat(X_max-X_min,[n_data,1]));

%% median correlation
W_1 = coeff1.W_1;
W_2 = coeff1.W_2;
W_out = coeff1.W_out;

b_1 = coeff1.b_1;
b_2 = coeff1.b_2;
b_out = coeff1.b_out;

Y1 = tanh(tanh(X_norm*W_1+b_1)*W_2+b_2)*W_out+b_out;
rho50 = reshape(Y1,n_row,n_col);

%% standard deviation (sigma_z)
W_1 = coeff2.W_1;
W_2 = coeff2.W_2;
W_out = coeff2.W_out;

b_1 = coeff2.b_1;
b_2 = coeff2.b_2;
b_out = coeff2.b_out;

Y2 = tanh(tanh(X_norm*W_1+b_1)*W_2+b_2)*W_out+b_out;
sigma_z = reshape(exp(Y2),n_row,n_col);

