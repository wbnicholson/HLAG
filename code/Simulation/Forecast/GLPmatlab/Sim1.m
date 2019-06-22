% Simulation script Scenario 1
clear;

% Load data
format long g;
load sim1.dat;

% Sim settings
Nsim = 500; % number of simulations
n = 100; % time series length
k = 45; % number of time series
T1=n; % first time point forecast evaluation window
T2=n; % last time point forecast evaluation window
p=5; % number of lags

% Store Forecast Results
MSFEs_GLP_sim1 = zeros(Nsim, 1); % Mean Squared Forecast Errors GLP


for i=1:Nsim
    
    
    Y=sim1(((n)*(i-1)+1):(i*n), 1:end);
    MSFE_GLP_sim1=zeros(T2-T1+1,1);

    for j=T1:T2
        
        Yraw=Y(1:j, 1:end);
        Ytest=Yraw(end,1:end);

        
        % GLP paper
        tic;
        YrawGLP=Yraw(1:(n-1), 1:end);
        res = bvarGLP(YrawGLP, p); 
        toc;
        YhatGLP = res.postmax.forecast(:,1:end);
        MSFE_GLP_sim1(j-T1+1,1)=sum((YhatGLP-Ytest).^2);
        
    end

    MSFEs_GLP_sim1(i,1)= mean(MSFE_GLP_sim1);

    
end

mean(MSFEs_GLP_sim1)/k

