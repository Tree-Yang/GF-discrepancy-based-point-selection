function [rps,asgn_prob] = GFDiscr_RPS(distr_para,n_rps)

%GFDiscr_RPS - Description
%GF-discrepancy based representative point selection strategy

%by J. S. Yang
%date: 2020-01-02

%====================================================================
    %INUPT:
    %----------------------------
    %distr_para: a structure including dimension, distribution type
    %            distribution parameters
    %n_rps: number of representative to be generated

    %OUTPUT:
    %----------------------------
    %rps: representative point set
    %asgn_prob: assign probabilities
%=====================================================================

    %basic parameter
    n_dim   = distr_para.dim;
    distr_type_a = distr_para.type;
    distr_data_a = distr_para.data;
    if length(distr_type_a) ~= n_dim
        error('Some error detected in "distr_para"');
    end
    n_skip   = 1e3;
    n_leap   = 1e2;
    n_mcs    = 1e7;
    view_opt = 'off';  % turn it on to debug

    %initial points: sobol set
    p_sobol  = sobolset(n_dim,'Skip',n_skip,'Leap',n_leap);
    p_sobol  = scramble(p_sobol,'MatousekAffineOwen');
    p_ini    = net(p_sobol, n_rps);

    %view the point set
    if strcmpi(view_opt, 'on')
        figure;
        plot(p_ini(:,1), p_ini(:,2), 'o', 'LineWidth', 1, 'MarkerSize',6, ...
        'MarkerEdgeColor',[31,78,121]/255, 'MarkerFaceColor',[157,195,230]/255);
        xlabel('\itX_{\rm1}'); ylabel('\itX_{\rm2}');
        set(gca, 'FontSize', 12, 'FontName', 'Arial');
    end


    %transformation via inverse cumulative distribution function to the physics space
    p_1 = zeros(size(p_ini));
    parfor ii = 1:1:n_dim
        distr_type = distr_type_a{ii};
        distr_data = distr_data_a{ii};
        if sum(strcmpi(distr_type, {'norm','normal'})) >=1
            p_1(:,ii) = norminv(p_ini(:,ii),distr_data(1),distr_data(2));
        elseif sum(strcmpi(distr_type, {'lognorm','lognormal'})) >=1
            p_1(:,ii) = logninv(p_ini(:,ii),distr_data(1),distr_data(2));
        elseif sum(strcmpi(distr_type, {'unif','uniform'})) >=1
            p_1(:,ii) = unifinv(p_ini(:,ii),distr_data(1),distr_data(2));
        else
            error('Wrong type for distribution!');
        end
    end
    %view the point set
    if strcmpi(view_opt, 'on')
        figure;
        plot(p_1(:,1), p_1(:,2), 'o', 'LineWidth', 1, 'MarkerSize',6, ...
        'MarkerEdgeColor',[197,90,17]/255, 'MarkerFaceColor',[248,203,173]/255);
        xlabel('\itX_{\rm1}'); ylabel('\itX_{\rm2}');
        set(gca, 'FontSize', 12, 'FontName', 'Arial');
        %GF-discrepancy
        asgn_prob = AsgnProb(p_1,distr_para,n_mcs);
        discr = GFDiscr(p_1, distr_para, asgn_prob);
        fprintf('The GF-discrepancy of the initial point set is %15.6f\n', discr);
    end

    %transformation via empirical cumulative distribution function
    p_20 = zeros(size(p_1));
    parfor jj = 1:1:n_rps
        p_1_tmp    = p_1;
        ind_tmp    = heaviside(repmat(p_1_tmp(jj,:), n_rps, 1)-p_1_tmp);
        p_20(jj,:) = sum(ind_tmp)/n_rps;
    end
    p_2 = zeros(size(p_20));
    parfor ii = 1:1:n_dim
        distr_type = distr_type_a{ii};
        distr_data = distr_data_a{ii};
        if sum(strcmpi(distr_type, {'norm','normal'})) >=1
            p_2(:,ii) = norminv(p_20(:,ii),distr_data(1),distr_data(2));
        elseif sum(strcmpi(distr_type, {'lognorm','lognormal'})) >=1
            p_2(:,ii) = logninv(p_20(:,ii),distr_data(1),distr_data(2));
        elseif sum(strcmpi(distr_type, {'unif','uniform'})) >=1
            p_2(:,ii) = unifinv(p_20(:,ii),distr_data(1),distr_data(2));
        else
            error('Wrong type for distribution!');
        end
    end
    %assigned probability
    asgn_prob = AsgnProb(p_2,distr_para,n_mcs);
    %view the point set
    if strcmpi(view_opt, 'on')
        figure;
        plot(p_2(:,1), p_2(:,2), 'o', 'LineWidth', 1, 'MarkerSize',6, ...
        'MarkerEdgeColor',[84,130,53]/255, 'MarkerFaceColor',[197,224,180]/255);
        xlabel('\itX_{\rm1}'); ylabel('\itX_{\rm2}');
        set(gca, 'FontSize', 12, 'FontName', 'Arial');
        %GF-discrepancy
        discr = GFDiscr(p_2, distr_para, asgn_prob);
        fprintf('The GF-discrepancy of the point set after the first rearrangement is %15.6f\n', discr);
    end

    %final transformation
    p_30 = zeros(size(p_2));
    parfor jj = 1:1:n_rps
        p_2_tmp      = p_2;
        ind_tmp     = heaviside(repmat(p_2_tmp(jj,:), n_rps, 1)-p_2_tmp);
        p_30(jj, :) = asgn_prob' * ind_tmp;
    end
    p_3 = zeros(size(p_30));
    parfor ii = 1:1:n_dim
        distr_type = distr_type_a{ii};
        distr_data = distr_data_a{ii};
        if sum(strcmpi(distr_type, {'norm','normal'})) >=1
            p_3(:,ii) = norminv(p_30(:,ii),distr_data(1),distr_data(2));
        elseif sum(strcmpi(distr_type, {'lognorm','lognormal'})) >=1
            p_3(:,ii) = logninv(p_30(:,ii),distr_data(1),distr_data(2));
        elseif sum(strcmpi(distr_type, {'unif','uniform'})) >=1
            p_3(:,ii) = unifinv(p_30(:,ii),distr_data(1),distr_data(2));
        else
            error('Wrong type for distribution!');
        end
    end
    %assigned probability
    asgn_prob = AsgnProb(p_3,distr_para,n_mcs);
    %view the point set
    if strcmpi(view_opt, 'on')
        figure;
        plot(p_3(:,1), p_3(:,2), 'o', 'LineWidth', 1, 'MarkerSize',6, ...
        'MarkerEdgeColor',[191,144,0]/255, 'MarkerFaceColor',[255,242,204]/255);
        xlabel('\itX_{\rm1}'); ylabel('\itX_{\rm2}');
        set(gca, 'FontSize', 12, 'FontName', 'Arial');
        %GF-discrepancy
        discr = GFDiscr(p_3, distr_para, asgn_prob);
        fprintf('The GF-discrepancy of the point set after the second rearrangement is %15.6f\n', discr);
    end
    rps = p_3;
end