function discr = GFDiscr(rps, distr_para, asgn_prob)

% GFDiscr - Description
% Evaluate the generalized F-discrepancy

% by J. S. Yang
% Date: 2020-01-03

%====================================================================
%INUPT:
%--------------------------------------------------------------------
    %rps: representative point set
    %distr_para: a structure including dimension, distribution type
    %            distribution parameters
    %asgn_prob: assign probabilities

%OUTPUT:
%--------------------------------------------------------------------
    %discr: generalized F-discrepancy
%=====================================================================

    %basic parameters
    n_dim = distr_para.dim;
    n_rps = size(rps,1);
    distr_type_a = distr_para.type;
    distr_data_a = distr_para.data;

    %empirical cumulative distribution function
    cdf_ed = zeros(size(rps));
    for jj = 1:1:n_rps
        rps0 = rps;
        % ind_tmp       = heaviside(repmat(rps0(jj,:), n_rps, 1)-rps0);
        ind_tmp       = (repmat(rps0(jj,:), n_rps, 1)-rps0) >= 0.0;
        cdf_ed(jj, :) = asgn_prob' * ind_tmp;
    end
    
    %theortical cumulative distribution function
    cdf_d  = zeros(size(rps));
    parfor ii = 1:1:n_dim
        distr_type = distr_type_a{ii};
        distr_data = distr_data_a{ii};
        if sum(strcmpi(distr_type, {'norm','normal'})) >=1
            cdf_d(:,ii) = normcdf(rps(:,ii),distr_data(1),distr_data(2));
        elseif sum(strcmpi(distr_type, {'lognorm','lognormal'})) >=1
            cdf_d(:,ii) = logncdf(rps(:,ii),distr_data(1),distr_data(2));
        elseif sum(strcmpi(distr_type, {'unif','uniform'})) >=1
            cdf_d(:,ii) = unifcdf(rps(:,ii),distr_data(1),distr_data(2));
        else
            error('Wrong type for distribution!');
        end
    end
    discr0 = max(abs(cdf_d-cdf_ed));
    discr  = max(discr0);
end