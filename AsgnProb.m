function asgn_prob = AsgnProb(rps, distr_para, n_mcs)
% AsgnProb - Description
% Evaluate assign probabilities via Monte Carlo simulation

%by J. S. Yang
%date: 2020-01-02

%==============================================================
% INPUT: 
%--------------------------------------
% rps: representative point set
% distr_para: structure including number of dimension, types of
%             distribution, parameters of distribution
% n_mcs: number of samples in MCS in order to evaluate the assign probabilities

% OUTPUT:
%--------------------------------------
% asgn_prob: the vector of assigned probabilities
%==============================================================

    n_dim   = distr_para.dim;
    distr_type_a = distr_para.type;
    distr_data_a = distr_para.data;
    pap = zeros(n_mcs,n_dim);

    %generate and normalize random samples and normalize representative points
    rps1 = zeros(size(rps));
    parfor ii = 1:1:n_dim
        distr_type = distr_type_a{ii};
        distr_data = distr_data_a{ii};
        if sum(strcmpi(distr_type, {'norm','normal'})) >=1              %normal distribution
            %generate samples
            pap(:,ii)  = normrnd(distr_data(1),distr_data(2),n_mcs,1);
            %mean value and standard devriation
            rv_mean = distr_data(1);
            rv_stdv = distr_data(2);
            %normalize samples
            pap(:,ii)  = (pap(:,ii)-rv_mean) / rv_stdv;
            %normalize point set
            rps1(:,ii) = (rps(:,ii)-rv_mean) / rv_stdv;
        elseif sum(strcmpi(distr_type, {'lognorm','lognormal'})) >=1    %log-normal distribution
            %generate samples
            pap(:,ii)  = lognrnd(distr_data(1),distr_data(2),n_mcs,1);
            %mean value and standard devriation
            rv_mean = exp(distr_data(1) + distr_data(2)^2/2);
            rv_stdv = rv_mean * (exp(distr_data(2)^2)-1);
            %normalize samples
            pap(:,ii)  = (pap(:,ii)-rv_mean) / rv_stdv;
            %normalize samples
            rps1(:,ii) = (rps(:,ii)-rv_mean) / rv_stdv;
        elseif sum(strcmpi(distr_type, {'unif','uniform'})) >=1         %uniform distribution
            %generate samples
            pap(:,ii)  = unifrnd(distr_data(1),distr_data(2),n_mcs,1);
            %mean value and standard devriation
            rv_mean = (distr_data(1) + distr_data(2)) / 2;
            rv_stdv = (distr_data(2) - distr_data(1))^2 / 12;
            %normalize samples
            pap(:,ii)  = (pap(:,ii)-rv_mean) / rv_stdv;
            %normalize samples
            rps1(:,ii) = (rps(:,ii)-rv_mean) / rv_stdv;
        else
            error('Wrong type for distribution!');
        end
    end

    %Voronoi cells
    %==============================================================
    % two methods are included here to evaluate the distance. The second 
    % one is much better than the first one on efficiency
    % tic;
    %-------------------------------------------------------------
    % Method I: parfor - slow
    % Date : 2020-01-02
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % [n_rps,~] = size(rps);
    % tmp = zeros(n_mcs,1);
    % parfor kk = 1:1:n_mcs
    %     dis = sqrt(sum((repmat(pap(kk,:),n_rps,1) - rps).^2,2));
    %     [~,tmp(kk)] = min(dis); 
    % end
    %-------------------------------------------------------------
    % Method II: knnsearch - fast
    % Date : 2020-01-03
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    tmp = knnsearch(rps1,pap);
    % toc;
    %==============================================================

    %==============================================================
    %Two methods are included here to count the number of points, the 
    % second method has higer effiency.
    % tic
    %-------------------------------------------------------------
    % Method I: parfor - slow
    % Date: 2020-01-02
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % asgn_prob = zeros(n_rps,1);
    % parfor ll = 1:1:n_rps
    %     asgn_prob(ll) = length(find(tmp == ll));
    % end
    %-------------------------------------------------------------
    % Method II: tabulate - fast
    % Date : 2020-01-03
    % Update: 2023-09-25
    % Sometimes, some representative points cannot be indexed by `knnsearch`.
    % In other words, the not indexed representative points are too far 
    %       from the Monte Carlo points, such that they are not nearest neighboors
    %       of any Monte Carlo point.
    % Then, the length of asgn_prob will not be equal to the number of representative points.
    % This problem is fixed by considering the first column of the output of `tabulate` function
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % number of representative points
    [n_rps,~] = size(rps);
    % initialize a vector for assigned probability
    asgn_prob = zeros(n_rps, 1);
    % count the frequency number of representative points
    tmp0 = tabulate(tmp);
    asgn_prob(int64(tmp0(:,1))) = tmp0(:,2);
    % toc
    %==============================================================

    asgn_prob = asgn_prob / n_mcs;
end