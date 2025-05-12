function [rho50, sigma_z] = Model3_xSI_IM_MSAS(X1,X2,X3,X4)

% Created by Mao-Xin Wang (dr.maoxin.wang@gmail.com or maoxinwang@ust.hk)
% January 2025
%
% Model 3: Correlations between ASI,SI,or DSI and other IMs (xSI-IM pairs)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%
%   X1   =  scalar or matrix of I_(MS,IM1) (set to 1, 2, ..., or 9 when MS IM denotes
%                 ASI, SI, DSI, PGA, PGV, AI, CAV, Ds5-75, or Ds5-95, respectively)
%   X2   =  scalar or matrix of I_(AS,IM2) (set to 1, 2, ..., or 9 when AS IM denotes
%                 ASI, SI, DSI, PGA, PGV, AI, CAV, Ds5-75, or Ds5-95, respectively)
%   X3   =  scalar or matrix of ksi (damping ratio, in percentage)
%   X4   =  scalar or matrix of I_HV (set to 1, 2, 3, or 4 for 
%           MS_H-AS_H, MS_H-AS_V, MS_V-AS_H,or MS_V-AS_V, respectively)
%          (Note: the above inputs must be in the same dimension)
%
% OUTPUT
%
%   rho50     =  scalar or matrix of median correlation
%   sigma_z   =  scalar or matrix of standard deviation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load coeff_Model3
X_min = [1 	1 	log(0.5) 	1];
X_max = [9 	9 	log(30) 	4];
[n_row,n_col] = size(X1);
n_data = n_row*n_col;
X1_temp = reshape(X1,n_data,1);
X2_temp = reshape(X2,n_data,1);
X3_temp = reshape(log(X3),n_data,1);
X4_temp = reshape(X4,n_data,1);
X_norm = ([X1_temp,X2_temp,X3_temp,X4_temp]-repmat(X_min,[n_data,1]))./(repmat(X_max-X_min,[n_data,1]));

%% median correlation
W_1 = coeff1.W_1;
W_out = coeff1.W_out;

b_1 = coeff1.b_1;
b_out = coeff1.b_out;

Y1 = tanh(X_norm*W_1+b_1)*W_out+b_out;
rho50 = reshape(Y1,n_row,n_col);

%% standard deviation (sigma_z)
W_1 = coeff2.W_1;
W_out = coeff2.W_out;

b_1 = coeff2.b_1;
b_out = coeff2.b_out;

Y2 = tanh(X_norm*W_1+b_1)*W_out+b_out;
sigma_z = reshape(exp(Y2),n_row,n_col);

