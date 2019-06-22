% VAR using the Gibbs sampler, based on independent Normal Wishar prior
%--------------------------------------------------------------------------
% Bayesian estimation, prediction and impulse response analysis in VAR
% models using the Gibbs sampler. 

clear all;
randn('seed',2); %#ok<RAND>
rand('seed',2); %#ok<RAND>
% Load data
format long g;
load Ysmallmedium_testh1.dat;
load Ysmallmedium_trainh1.dat;

check = size(Ysmallmedium_trainh1);
rolsize = check(1)/76;

%----------------------------PRELIMINARIES---------------------------------
p=4; % number of lags
k=check(2); % number of time series
% Define specification of the VAR model
constant = 1;        % 1: if you desire intercepts, 0: otherwise 
forecasting = 1;     % 1: Compute h-step ahead predictions, 0: no prediction
forecast_method = 1; % 0: Direct forecasts 
                     % 1: Iterated forecasts
repfor = 10;         % Number of times to obtain a draw from the predictive 
                     % density, for each generated draw of the parameters                     
h = 1;               % Number of forecast periods



% Gibbs-related preliminaries
nsave = 1000;         % Final number of draws to save
nburn = 0.1*nsave;         % Draws to discard (burn-in)
ntot = nsave + nburn;  % Total number of draws
it_print = 1000;      % Print on the screen every "it_print"-th iteration

%========================== Prior Stuff ================================
%==========================================================================

%% -----------------Prior hyperparameters for bvar model
theta=[0.1 0.5 100 2];       % hyperparameters of Minnesota prior: 
                             % [lambda1 lambda2 int lambda3], int is the 
                             % prior on the intercept. lambda1, lambda2
                             % and lambda3 are as in equation (42) with
                             % lambda1 the overall shrinkage, lambda2 the
                             % cross srhinkage and lambda 3 the lag decay
                             % (quadratic if =2). Note lambda2~=1 implies
                             % the prior becomes asymmetric across eqation,
                             % so this would not be implementable in the
                             % standard conjugate setup. 
Minn_pmean = 0;              % Prior mean of the 1-st own lag for each 
                             % equation. For nonstationary variables, this
                             % is usually set to 1. For transformed
                             % stationary variables this is set to 0.  
                             

MSFEs_CCM_smallmedium_h1 = zeros(76, k);

for isim=1:76
    
disp(sprintf('now running isim= %d', isim));
        
Yraw = Ysmallmedium_trainh1(((isim-1)*rolsize+1):(rolsize*isim), 1:end);

if h>1
    for ih=1:(h-1)
        Yraw = [Yraw(1:end, 1:end) ; zeros(1, k)]; 
    end
    Yraw = [Yraw(1:end, 1:end) ; Ysmallmedium_testh1(isim,1:end)]; 
end

if h==1
  Yraw = [Yraw(1:end, 1:end) ; Ysmallmedium_testh1(isim,1:end)];   
end


%--------------------------DATA HANDLING-----------------------------------
% Get initial dimensions of dependent variable
[Traw M] = size(Yraw);

% The model specification is different when implementing direct forecasts,
% compared to the specification when computing iterated forecasts.
if forecasting==1
    if h<=0
        error('You have set forecasting, but the forecast horizon choice is wrong')
    end    

    % Now create VAR specification according to forecast method
    if forecast_method==0       % Direct forecasts
        Y1 = Yraw(h+1:end,:);
        Y2 = Yraw(2:end-h,:);
        Traw = Traw - h - 1;
    elseif forecast_method==1   % Iterated forecasts
        Y1 = Yraw;
        Y2 = Yraw;
    else
        error('Wrong choice of forecast_method')
    end
else
   Y1 = Yraw;
   Y2 = Yraw;
end

% Generate lagged Y matrix. This will be part of the X matrix
Ylag = mlag2(Y2,p); % Y is [T x M]. ylag is [T x (Mp)]

% Now define matrix X which has all the R.H.S. variables (constant, lags of
% the dependent variable and exogenous regressors/dummies)
if constant
    X1 = [ones(Traw-p,1) Ylag(p+1:Traw,:)];
else
    X1 = Ylag(p+1:Traw,:);  %#ok<UNRCH>
end

% Get size of final matrix X
[Traw3 K] = size(X1);


% Form Y matrix accordingly
% Delete first "LAGS" rows to match the dimensions of X matrix
Y1 = Y1(p+1:Traw,:); % This is the final Y matrix used for the VAR

% Traw was the dimesnion of the initial data. T is the number of actual 
% time series observations of Y and X (we lose the p-lags)
T = Traw - p;

%========= FORECASTING SET-UP:
% Now keep also the last "h" or 1 observations to evaluate (pseudo-)forecasts
if forecasting==1
    Y_pred = zeros(nsave*repfor,M); % Matrix to save prediction draws
    PL =zeros(nsave,1);             % Matrix to save Predictive Likelihood
    
    if forecast_method==0  % Direct forecasts, we only need to keep the 
        Y = Y1(1:end-1,:);                             % last observation
        X = X1(1:end-1,:);
        T = T - 1;
    else              % Iterated forecasts, we keep the last h observations
        Y = Y1(1:end-h,:);
        X = X1(1:end-h,:);
        T = T - h;
    end
else
    Y = Y1;
    X = X1;
%     Z = Z1;
end


%-----------------------------PRELIMINARIES--------------------------------
SIGMA = eye(M); 

% Storage space for posterior draws
ALPHA_draws = zeros(nsave,K,M);
SIGMA_draws = zeros(nsave,M,M);

% Prior on conditional mean coefficients, use Minnesota setup
ARresid=[];
for i=1:M
    yt_0=[ones(T-1,1) Y(1:end-1,i)];
    yt_1=Y(2:end,i);
    ARresid(:,i)=yt_1-yt_0*(yt_0\yt_1); %#ok<SAGROW>
end
AR_s2= diag(diag(ARresid'*ARresid))./(T-2);

Pi_pm=zeros(M*(K-1),1); Pi_pv=ones(M*(K-1), 1); co=0;
for i=1:M
    sigma_const(i)=AR_s2(i,i)*theta(3); %#ok<SAGROW> % this sets the prior variance on the intercept 
    for l=1:p;
        for j=1:M
            co=co+1;
            if (i==j)
                Pi_pv(co,1)=theta(1)/(l^theta(4)); % prior variance, own lags
                if l==1; Pi_pm(co)=Minn_pmean; end; % this sets the prior means for the first own lag coefficients.
            else
                Pi_pv(co,1)=(AR_s2(i,i)/AR_s2(j,j)*theta(1)*theta(2)/(l^theta(4))); % prior variance, cross-lags
            end
        end
    end
end

OMEGA_pai_vec = vec([sigma_const;reshape(Pi_pv,K-1, M)]);
MU_pai      = [zeros(1,M);reshape(Pi_pm,K-1,M)];                   % prior mean of Pai=ALPHA 

PAI = zeros(K, M);                                         % pre-allocate space for PAI
comp=[eye(M*(p-1)),zeros(M*(p-1),M)];                     % companion form
iV_vec = 1./OMEGA_pai_vec;
iVb_prior=iV_vec.*vec(MU_pai);    % inverses of prior matrices

 % Hyperparameters on inv(SIGMA) ~ W(v_prior,inv(S_prior))
 v_prior = M;
 S_prior = eye(M);
 inv_S_prior = inv(S_prior);  
    

%========================== Start Sampling ================================
%==========================================================================
tic;
for irep = 1:ntot  %Start the Gibbs "loop"

    % This is the only new step (triangular algorithm).
    [L,D] = ldl(SIGMA);
	invA_ = L;
    sqrt_ht = repmat(sqrt(diag(D)'), T, 1);
    PAI = triang(Y, X, M, K, T, invA_, sqrt_ht, iV_vec, iVb_prior);  
    ALPHA = PAI;

    % Posterior of SIGMA|ALPHA,Data ~ iW(inv(S_post),v_post)
    v_post = T + v_prior;
    S_post = S_prior + (Y - X*ALPHA)'*(Y - X*ALPHA);
    SIGMA = inv(wish(inv(S_post),v_post));% Draw SIGMA

    % Store results  
    if irep > nburn               
        %=========FORECASTING:
        if forecasting==1
            if forecast_method == 0   % Direct forecasts
                Y_temp = zeros(repfor,M);
                % compute 'repfor' predictions for each draw of ALPHA and SIGMA
                for ii = 1:repfor
                    X_fore = [1 Y(T,:) X(T,2:M*(p-1)+1)];
                    % Forecast of T+1 conditional on data at time T
                    Y_temp(ii,:) = X_fore*ALPHA + randn(1,M)*chol(SIGMA);
                end
                % Matrix of predictions
                Y_pred(((irep-nburn)-1)*repfor+1:(irep-nburn)*repfor,:) = Y_temp;
                % Predictive likelihood
                PL(irep-nburn,:) = mvnpdf(Y1(T+1,:),X(T,:)*ALPHA,SIGMA);
                if PL(irep-nburn,:) == 0
                    PL(irep-nburn,:) = 1;
                end
            elseif forecast_method == 1   % Iterated forecasts
                % The usual way is to write the VAR(p) model in companion
                % form, i.e. as VAR(1) model in order to estimate the
                % h-step ahead forecasts directly (this is similar to the 
                % code we use below to obtain impulse responses). Here we 
                % just iterate over h periods, obtaining the forecasts at 
                % T+1, T+2, ..., T+h iteratively.
                
                for ii = 1:repfor
                    % Forecast of T+1 conditional on data at time T
                    X_fore = [1 Y(T,:) X(T,2:M*(p-1)+1)];
                    Y_hat = X_fore*ALPHA + randn(1,M)*chol(SIGMA);
                    Y_temp = Y_hat;
                    X_new_temp = X_fore;
                    for i = 1:h-1  % Predict T+2, T+3 until T+h                   
                        if i <= p
                            % Create matrix of dependent variables for                       
                            % predictions. Dependent on the horizon, use the previous                       
                            % forecasts to create the new right-hand side variables
                            % which is used to evaluate the next forecast.                       
                            X_new_temp = [1 Y_hat X_fore(:,2:M*(p-i)+1)];
                            % This gives the forecast T+i for i=1,..,p                       
                            Y_temp = X_new_temp*ALPHA + randn(1,M)*chol(SIGMA);                       
                            Y_hat = [Y_hat Y_temp];
                        else
                            X_new_temp = [1 Y_hat(:,1:M*p)];
                            Y_temp = X_new_temp*ALPHA + randn(1,M)*chol(SIGMA);
                            Y_hat = [Y_hat Y_temp];
                        end
                    end %  the last value of 'Y_temp' is the prediction T+h
                    Y_temp2(ii,:) = Y_temp;
                end
                % Matrix of predictions               
                Y_pred(((irep-nburn)-1)*repfor+1:(irep-nburn)*repfor,:) = Y_temp2;
                % Predictive likelihood
                PL(irep-nburn,:) = mvnpdf(Y1(T+h,:),X_new_temp*ALPHA,SIGMA);
                if PL(irep-nburn,:) == 0
                    PL(irep-nburn,:) = 1;
                end
            end
        end % end forecasting
        %=========Forecasting ends here
        
    
               
        %----- Save draws of the parameters
        ALPHA_draws(irep-nburn,:,:) = ALPHA;
        SIGMA_draws(irep-nburn,:,:) = SIGMA;

    end % end saving results
       
end %end the main Gibbs for loop
%====================== End Sampling Posteriors ===========================

% mean prediction and log predictive likelihood
if forecasting == 1
    Y_pred_mean = mean(Y_pred,1);

%This are the true value of Y at T+h:
if forecast_method == 0
    true_value = Y1(T+1,:);
elseif forecast_method == 1
    true_value = Y1(T+h,:);
end
MSFE = (Y_pred_mean - true_value).^2;
end

clc;
toc;

MSFEs_CCM_smallmedium_h1(isim,1:end) = MSFE;

end
