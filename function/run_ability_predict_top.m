function [coef,coef0]= run_ability_predict_top(n,x,y)

    xTrain = x;
    yTrain = y;

    % ==== Set some parameters for prediction ==== %
    %loop over the two independent samples

    coef0=zeros(1,1);
    lambda=zeros(1,1);
    penalty='ridge';%'lasso','ridge''kridge';
    %Range for lambda
    L='default';
    %Max number lambda values evaluated
    J = 100;
    %Number of cross-validations
    ncvals = 20;
    %Search strategy {'gridsearch','bayes'}
    opt_method='gridsearch';
    %Show a plot of the error as a function of lambda
    ShowPlot=0;

    if strcmp(penalty,'kridge')
        nSubjects=size(xTrain,1);
        coef=zeros(nSubjects,1); %beta coeffients
    else
        coef=zeros(size(xTrain,2),1); % number of edges x 2 groups Coefficients of edges
    end

    % Train is used to train the model and test is used to test
    x=xTrain; y=yTrain;

    % Fit model
    if strcmp(penalty,'lasso')
        solver={'sparsa'};
    elseif strcmp(penalty,'ridge')
        solver={'asgd','lbfgs'};
    elseif strcmp(penalty,'kridge')
        solver={'asgd','lbfgs'};
        xtest=corr(xtest',x'); %subject x subject %Reference back to kernel
        x=corr(x');%subject x subject
        penalty='ridge';
    end

    params=hyperparameters('fitrlinear',x,y);

    if ~strcmp(L,'default')
        params(1).Range=L;
    end

        params(2).Optimize=false; params(3).Optimize=false;
        % cognitive measures
        [mdl,~,hyper]=fitrlinear(x',y,... % transpose x->x':significant reduction in optimization-execution time
            'Learner','leastsquares',...
            'Regularization',penalty,...
            'OptimizeHyperparameters',params,...
            'Verbose',0,...
            'Solver',solver,...  %default is 'sgd'
            'ObservationsIn','columns',... % transpose x->x':significant reduction in optimization-execution time
            'PostFitBias',true,...         %default is false
            'PassLimit',10,...             %default is 1 (increasing will increase run time)
            'HyperparameterOptimizationOptions',struct('Kfold',ncvals,...
            'Optimizer',opt_method,...
            'NumGridDivisions',J,...
            'ShowPlots',false,...
            'Repartition',false,...
            'MaxObjectiveEvaluations',J));

        coef(:,1)=mdl.Beta; coef0(1)=mdl.Bias; lambda(1)=mdl.Lambda;
end