function PAI=triang(Y,X,N,K,T,invA_,sqrt_ht,iV_vec,iVb_prior)
% =========================================================================
% Performs a draw from the conditional posterior of the VAR conditional
% mean coefficients by using the triangular algorithm. The
% triangularization achieves computation gains of order N^2 where N is the
% number of variables in the VAR. Carriero, Clark and Marcellino (2015),
% Large Vector Autoregressions with stochastic volatility and flexible
% priors. 
%
% The model is:
%
%     Y(t) = Pai(L)Y(t-1) + v(t); Y(t) is Nx1; t=1,...,T, L=1,...,p.
%     v(t) = inv(A)*(LAMBDA(t)^0.5)*e(t); e(t) ~ N(0,I);
%                _                                         _
%               |    1          0       0       ...      0  |
%               |  A(2,1)       1       0       ...      0  |
%      A =      |  A(3,1)     A(3,2)    1       ...      0  |
%               |   ...        ...     ...      ...     ... |
%               |_ A(N,1)      ...     ...   A(N,N-1)    1 _|
%
%     Lambda(t)^0.5 = diag[sqrt_h(1,t)  , .... , sqrt_h(N,t)];
%
% INPUTS
% Data and pointers:
% Y     = (TxN) matrix of data appearing on the LHS of the VAR
% X     = (TxK) matrix of data appearing on the RHS of the VAR
% N     = scalar, #of variables in VAR 
% K     = scalar, #of regressors (=N*p+1)  
% T     = scalar, #of observations
% The matrix X needs to be ordered as: [1, y(t-1), y(t-2),..., y(t-p)]
% 
% Error variance stuff:
% invA_   = (NxN) inverse of lower triangular covariance matrix A
% sqrt_ht = (TxN) time series of diagonal elements of volatility matrix 
% For a homosckedastic system, with Sigma the error variance, one can
% perform the LDL decomposition (command [L,D]=LDL(Sigma)) and set inv_A=L
% and sqrt_ht=repmat(sqrt(diag(D)'),T,1). 
%
% Priors:
% iV          = (NKxNK) precision matrix for VAR coefficients 
% iVB_prior   = (NKx1) (prior precision)*(prior mean)
% Note 1:iV is the inverse of the prior matrix and iVB_prior is the product
% of iV and the prior mean vector, which both need to be computed only once,
% before the start of the main MCMC loop.
% Note 2:in this code, iV is assumed block-diagonal. This implies that the
% prior is independent across equations. This includes most of the priors 
% usually considered, including the Minnesota one.  To use a non-block
% diagonal iV one needs to modify the code using the recursions illustrated 
% in equations (37) and (38).  
%
% OUTPUT
% One draw from (PAI|A,Lambda,data)
% PAI=[Pai(0), Pai(1), ..., Pai(p)].
% =========================================================================
PAI=zeros(K,N);
for j=1:N
    
    % iteratively compute previous equation residuals and rescale the data.
    % *** This is the sub-loop giving O(N^2) computational gains *** 
    Yr= Y(:,j);  
    if j==1; resid=zeros(T,N);
    else for l=1:j-1; Yr= Yr-invA_(j,l)*(sqrt_ht(:,l).*resid(:,l)); end
    end
    Yr= Yr./sqrt_ht(:,j); Xr= X./repmat(sqrt_ht(:,j),1,K);
    
    % index to select equation j coefficients from general specification
    index=K*(j-1)+1:(K*(j-1)+K);
    
    % perform the draw of PAI(equation j)|Pai(equations 1,...,j-1) 
    V_post = (diag(iV_vec(index)) + Xr'*Xr)\eye(K); % equation 33
    b_post = V_post*(iVb_prior(index) + Xr'*Yr); % equation 32
    pai_j  = b_post + chol(V_post,'lower')*randn(K,1); % draw
    
    % The 2 lines below perform a draw from pai_j using a formula similar
    % to that in footnote 5, but equation by equation. These 2 lines are
    % based on inv(Chol_precision)*inv(Chol_precision)'=chol(V_post,'lower')*chol(V_post,'lower')'
    % Using these 2 lines will give further (but minor) speed improvements. 
    % Chol_precision = chol(iV(index,index) + Xr'*Xr);
    % pai_j = Chol_precision\((Chol_precision'\(iVb_prior(index) + Xr'*Yr)) + randn(K,1));
    
    % save residuals and present draw of pai_j
    resid(:,j)=Yr-Xr*pai_j; % store residuals to be removed in next equation
    PAI(:,j) = pai_j;       % store the draw of PAI  
    
end
            
