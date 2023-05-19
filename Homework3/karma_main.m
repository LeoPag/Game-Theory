clear;
close all;
clc;

%% Parameters
% Urgency
param.U = [1; 10];
param.n_u = length(param.U);
param.prob_up_u = [0.6; 0.4];

% Karma
param.k_max = 20;
param.K = (0 : param.k_max).';
param.n_k = length(param.K);
param.k_bar = 10;
param.i_k_bar = find(param.K == param.k_bar);

% Outcomes
param.O = [0; 1];   % 0 := ego agent wins, 1 := ego agent yields
param.n_o = length(param.O);

% Algorithm parameter
param.alpha = 0.8;
param.pi_tol = 1e-6;
param.d_tol = 1e-6;
param.V_tol = 1e-10;
param.V_max_iter = 10000;
param.max_iter = 2000;
param.lambda = 1000;
param.dt = 0.05 : (param.pi_tol - 0.05) / (param.max_iter - 1) : param.pi_tol;

%% Pre-computations for speed
% Expected stage cost
cost_down_u_o = param.U * param.O.';

% Probability of winning/yielding given bids
prob_down_b_bj_up_o = zeros(param.n_k, param.n_k, param.n_o);
for i_b = 1 : param.n_k
    for i_bj = 1 : param.n_k
        if i_b > i_bj
            prob_down_b_bj_up_o(i_b,i_bj,1) = 1;
        elseif i_b < i_bj
            prob_down_b_bj_up_o(i_b,i_bj,1) = 0;
        else
            prob_down_b_bj_up_o(i_b,i_bj,1) = 0.5;
        end
        prob_down_b_bj_up_o(i_b,i_bj,2) = 1 - prob_down_b_bj_up_o(i_b,i_bj,1);
    end
end

% Probability of next karma given bids
prob_down_k_b_bj_up_kn_pre_tax = zeros(param.n_k, param.n_k, param.n_k, param.n_k);
for i_k = 1 : param.n_k
    for i_b = 1 : i_k
        b = param.K(i_b);
        for i_bj = 1 : param.n_k
            bj = param.K(i_bj);

            if i_b > i_bj   % Ego agent wins
                p_win = get_payments(b, bj);
                assert(p_win <= b, 'Error: payment is larger than bid');
                assert(mod(p_win, 1) == 0, 'Error: payment is non-integer');
                i_kn = min([i_k - p_win, param.n_k]);
                prob_down_k_b_bj_up_kn_pre_tax(i_k,i_b,i_bj,i_kn) = 1;

            elseif i_b < i_bj   % Ego agent yields
                [~, p_yield] = get_payments(bj, b);
                assert(p_yield <= b, 'Error: payment is larger than bid');
                assert(mod(p_yield, 1) == 0, 'Error: payment is non-integer');
                i_kn = min([i_k - p_yield, param.n_k]);
                prob_down_k_b_bj_up_kn_pre_tax(i_k,i_b,i_bj,i_kn) = 1;

            else    % Tie
                [p_win, p_yield] = get_payments(b, bj);
                assert(p_win <= b && p_yield <= bj, 'Error: payment is larger than bid');
                assert(mod(p_win, 1) == 0 && mod(p_yield, 1) == 0, 'Error: payment is non-integer');

                % Win with probability 0.5
                i_kn_win = min([i_k - p_win, param.n_k]);
                prob_down_k_b_bj_up_kn_pre_tax(i_k,i_b,i_bj,i_kn_win) = 0.5;

                % Yield with probability 0.5
                i_kn_yield = min([i_k - p_yield, param.n_k]);
                prob_down_k_b_bj_up_kn_pre_tax(i_k,i_b,i_bj,i_kn_yield) = prob_down_k_b_bj_up_kn_pre_tax(i_k,i_b,i_bj,i_kn_yield) + 0.5;
            end
        end
    end
end

% Surplus payment given bids (to be redistributed)
p_surplus_down_b_win_b_yield = zeros(param.n_k, param.n_k);
for i_b_win = 1 : param.n_k
    b_win = param.K(i_b_win);
    for i_b_yield = 1 : i_b_win
        b_yield = param.K(i_b_yield);
        [p_win, p_yield] = get_payments(b_win, b_yield);
        assert(p_win + p_yield >= 0, 'Error: surplus payment is negative');
        p_surplus_down_b_win_b_yield(i_b_win,i_b_yield) = p_win + p_yield;
    end
end

% Tax
tax_down_k = zeros(param.n_k, 1);
prob_down_k_pre_tax_up_k_post_tax = zeros(param.n_k, param.n_k);
for i_k = 1 : param.n_k
    tax = get_karma_tax(param.K(i_k));
    assert(tax <= param.K(i_k) && tax >= 0, 'Invalid karma tax');
    
    tax_low = floor(tax);
    i_k_post_tax_low = i_k - tax_low;
    if tax_low == tax
        prob_down_k_pre_tax_up_k_post_tax(i_k,i_k_post_tax_low) = 1;
    else
        tax_high = tax_low + 1;
        i_k_post_tax_high = i_k - tax_high;
        prob_down_k_pre_tax_up_k_post_tax(i_k,i_k_post_tax_low) = tax_high - tax;
        prob_down_k_pre_tax_up_k_post_tax(i_k,i_k_post_tax_high) = tax - tax_low;
    end
    tax_down_k(i_k) = tax;
end

%% Initialization
% Initial policy
pi_down_u_k_up_b = zeros(param.n_u, param.n_k, param.n_k);
for i_u = 1 : param.n_u
    for i_k = 1 : param.n_k
        % Bid 0.5 * u / u_max * k (~ 'bid half if urgent')
        b = round(0.5 * param.U(i_u) / param.U(end) * param.K(i_k));
        i_b = param.K == b;
        pi_down_u_k_up_b(i_u,i_k,i_b) = 1;
    end
end

% Initial state distribution
% Initially all agents have average karma k_bar
d_up_k = zeros(param.n_k, 1);
d_up_k(param.i_k_bar) = 1;
d_up_u_k = param.prob_up_u * d_up_k.';

%% Memory allocations for main algorithm
iter = 1;
pi_error = inf;
d_error = inf;
V_down_u_k = zeros(param.n_u, param.n_k);
V_down_u_k_next = zeros(param.n_u, param.n_k);
Q_down_u_k_b = zeros(param.n_u, param.n_k, param.n_k);
br_pi_down_u_k_up_b = zeros(param.n_u, param.n_k, param.n_k);

% Display status
fprintf('%% COMPUTING NASH EQULIBRIUM, PLEASE WAIT %%\n\n');
pause(1);

%% Main algorithm loop
while (pi_error > param.pi_tol || d_error > param.d_tol) && iter <= param.max_iter
    % Distribution of others' karma + bids
    prob_up_kj_bj = einsum('ijk,ij->jk', pi_down_u_k_up_b, d_up_u_k);
    
    % Distribution of others' bids
    prob_up_bj = sum(prob_up_kj_bj, 1).';

    % Probability of winning given bid
    prob_down_b_up_o = einsum('ijk,jl->ikl', prob_down_b_bj_up_o, prob_up_bj);

    % Stage rewards
    reward_down_u_b = -einsum('ij,kj->ik', cost_down_u_o, prob_down_b_up_o);
    R_down_u_k = einsum('ijk,ik->ij', pi_down_u_k_up_b, reward_down_u_b);
    
    % Probability of next karma given bid before tax and redistribution
    prob_down_k_b_up_kn_pre_tax = einsum('ijkl,km->ijl', prob_down_k_b_bj_up_kn_pre_tax, prob_up_bj);
    
    % Average surplus
    prob_up_kn_pre_tax = einsum('ij,ijk->k', prob_up_kj_bj, prob_down_k_b_up_kn_pre_tax);
    p_surplus_bar = param.k_bar - dot(prob_up_kn_pre_tax, param.K);
    
    % Tax
    tax_bar = dot(prob_up_kn_pre_tax, tax_down_k);
    prob_down_k_b_up_kn_pre_red = einsum('ijk,kl->ijl', prob_down_k_b_up_kn_pre_tax, prob_down_k_pre_tax_up_k_post_tax);
    
    % Redistribution
    redist = p_surplus_bar + tax_bar;
    redist_low = floor(redist);
    redist_high = redist_low + 1;
    prob_up_kn_pre_red = einsum('ij,ijk->k', prob_up_kj_bj, prob_down_k_b_up_kn_pre_red);
    saturation_adjustment_term = 0;
    for i = 1 : redist_low
        saturation_adjustment_term = saturation_adjustment_term + prob_up_kn_pre_red(param.n_k - i) * i;
    end
    redist_adjusted = (redist - saturation_adjustment_term) / sum(prob_up_kn_pre_red(1:param.n_k-redist_high));
    redist_adjusted_low = floor(redist_adjusted);
    redist_adjusted_high = redist_adjusted_low + 1;
    prob_down_kn_pre_red_up_kn = zeros(param.n_k, param.n_k);
    for i_k = 1 : param.n_k
        i_kn_low = min([i_k + redist_adjusted_low, param.n_k]);
        i_kn_high = min([i_k + redist_adjusted_high, param.n_k]);
        prob_down_kn_pre_red_up_kn(i_k,i_kn_low) = redist_adjusted_high - redist_adjusted;
        prob_down_kn_pre_red_up_kn(i_k,i_kn_high) = prob_down_kn_pre_red_up_kn(i_k,i_kn_high) + redist_adjusted - redist_adjusted_low;
    end
    prob_down_k_b_up_kn = einsum('ijk,kl->ijl', prob_down_k_b_up_kn_pre_red, prob_down_kn_pre_red_up_kn);
    
    % Stochastic matrix
    prob_down_k_b_up_un_kn = einsum('ijk,lm->ijlk', prob_down_k_b_up_kn, param.prob_up_u);
    prob_down_u_k_b_up_un_kn = einsum('ijkl,mn->mijkl', prob_down_k_b_up_un_kn, ones(param.n_u, 1));
    P_down_u_k_up_un_kn = einsum('ijk,ijklm->ijlm', pi_down_u_k_up_b, prob_down_u_k_b_up_un_kn);

    % Value function
    V_iter = 1;
    V_error = inf;
    while V_error > param.V_tol && V_iter <= param.V_max_iter
        if param.alpha == 1   % alpha = 1 is treated as an average cost per stage problem
            rel_V = R_down_u_k(1,param.i_k_bar) + einsum('ijkl,kl->ij', P_down_u_k_up_un_kn(1,param.i_k_bar,:,:), V_down_u_k);
            V_down_u_k_next = R_down_u_k + einsum('ijkl,kl->ij', P_down_u_k_up_un_kn, V_down_u_k) - rel_V;
        else
            V_down_u_k_next = R_down_u_k + param.alpha * einsum('ijkl,kl->ij', P_down_u_k_up_un_kn, V_down_u_k);
        end
        V_error = norm(reshape(V_down_u_k - V_down_u_k_next, 1, []), inf);
        V_down_u_k = V_down_u_k_next;
        V_iter = V_iter + 1;
    end

    % Single-stage deviation rewards
    r_down_u_k_b = einsum('ij,kl->ikj', reward_down_u_b, ones(param.n_k, 1));
    future_r_down_u_k_b = einsum('ijklm,lm->ijk', prob_down_u_k_b_up_un_kn, V_down_u_k);
    Q_down_u_k_b = r_down_u_k_b + param.alpha * future_r_down_u_k_b;

    % Perturbed best response policy
    lambda_Q_down_u_k_b = param.lambda * Q_down_u_k_b;
    for i_k = 1 : param.n_k
        br_pi_down_u_k_up_b(:,i_k,1:i_k) = exp(lambda_Q_down_u_k_b(:,i_k,1:i_k) - max(lambda_Q_down_u_k_b(:,i_k,1:i_k), [], 3)); % Subtract max for numerical stability
    end
    br_pi_down_u_k_up_b = br_pi_down_u_k_up_b ./ sum(br_pi_down_u_k_up_b, 3);

    % Next policy
    pi_down_u_k_up_b_next = (1 - param.dt(iter)) * pi_down_u_k_up_b + param.dt(iter) * br_pi_down_u_k_up_b;
    pi_error = norm(reshape(pi_down_u_k_up_b_next - pi_down_u_k_up_b, 1, []), inf);

    % Next state distribution
    d_up_u_k_next = einsum('ij,ijkl->kl', d_up_u_k, P_down_u_k_up_un_kn);
    d_up_u_k_next = d_up_u_k_next / sum(d_up_u_k_next(:));
    d_up_u_k_next = (1 - param.dt(iter)) * d_up_u_k + param.dt(iter) * d_up_u_k_next;
    d_error = norm(reshape(d_up_u_k_next - d_up_u_k, 1, []), inf);

    % Update policy & distribution candidates
    pi_down_u_k_up_b = pi_down_u_k_up_b_next;
    d_up_u_k = d_up_u_k_next;

    % Display status
    fprintf('Iteration %d policy error %f distribution error %f\n', iter, pi_error, d_error);

    % Increment iteration count
    iter = iter + 1;
end

%% Plot
screensize = get(groot, 'ScreenSize');
screenwidth = screensize(3);
screenheight = screensize(4);
load('RedColormap.mat');
figure(1);
fig = gcf;
fig.Position = [0, 0, 0.9 * screenwidth, 0.9 * screenheight];
n_subplots = 2 * param.n_u + 4;
n_rows = 2;
n_cols = ceil(n_subplots / n_rows);

% Policy plot
for i_u = 1 : param.n_u
    pi_mat = squeeze(br_pi_down_u_k_up_b(i_u,:,:));
    for i_k = 1 : param.n_k - 1
        pi_mat(i_k,i_k+1:param.n_k) = nan;
    end
    i_subplot = 2 * i_u - 1;
    subplot(n_rows, n_cols, [i_subplot, i_subplot + 1]);
    pi_plot = heatmap(param.K, param.K, pi_mat.', 'ColorbarVisible','off');
    pi_plot.YDisplayData = flipud(pi_plot.YDisplayData);
    pi_plot.Title = ['NE policy for u = ', num2str(param.U(i_u))];
    pi_plot.XLabel = 'Karma';
    pi_plot.YLabel = 'Bid';
    pi_plot.FontName = 'Ubuntu';
    pi_plot.FontSize = 20;
    if exist('RedColormap', 'var')
        pi_plot.Colormap = RedColormap;
    end
    pi_plot.ColorLimits = [0 1];
    pi_plot.CellLabelColor = 'none';
    pi_plot.GridVisible = false;
    for i_k = 1 : param.n_k
        if mod(i_k - 1, 5) ~= 0
            pi_plot.XDisplayLabels{i_k} = '';
            pi_plot.YDisplayLabels{i_k} = '';
        end
    end
end

% Stationary karma distribution plot
i_subplot = i_subplot + 2;
subplot(n_rows, n_cols, [i_subplot, i_subplot + 1]);
bar(param.K, sum(d_up_u_k, 1));
axis tight;
axes = gca;
axes.Title.FontName = 'ubuntu';
axes.Title.String = 'NE karma distribution';
axes.Title.FontSize = 20;
axes.XAxis.FontSize = 20;
axes.XAxis.FontName = 'ubuntu';
axes.YAxis.FontSize = 20;
axes.YAxis.FontName = 'ubuntu';
axes.XLabel.FontName = 'ubuntu';
axes.XLabel.String = 'Karma';
axes.XLabel.FontSize = 20;
axes.YLabel.FontName = 'ubuntu';
axes.YLabel.String = 'Distribution';
axes.YLabel.FontSize = 20;

% Stationary bid distribution plot
i_subplot = i_subplot + 2;
subplot(n_rows, n_cols, i_subplot);
bar(param.K, prob_up_bj, 'Facecolor', '#D95319');
axis tight;
axes = gca;
axes.Title.FontName = 'ubuntu';
axes.Title.String = 'NE bid distribution';
axes.Title.FontSize = 20;
axes.XAxis.FontSize = 20;
axes.XAxis.FontName = 'ubuntu';
axes.YAxis.FontSize = 20;
axes.YAxis.FontName = 'ubuntu';
axes.XLabel.FontName = 'ubuntu';
axes.XLabel.String = 'Bid';
axes.XLabel.FontSize = 20;
axes.YLabel.FontName = 'ubuntu';
axes.YLabel.String = 'Distribution';
axes.YLabel.FontSize = 20;

% Social cost
social_cost = -einsum('ijk,ijk->k', d_up_u_k, R_down_u_k);
optimal_social_cost = 0;
for i_u = 1 : param.n_u
    optimal_social_cost = optimal_social_cost + 0.5 * param.prob_up_u(i_u)^2 * param.U(i_u);
    for i_uj = i_u + 1 : param.n_u
        optimal_social_cost = optimal_social_cost + param.prob_up_u(i_u) * param.prob_up_u(i_uj) * param.U(i_u);
    end
end
social_costs = [social_cost, optimal_social_cost];
i_subplot = i_subplot + 1;
subplot(n_rows, n_cols, i_subplot);
bar(social_costs, 'Facecolor', '#EDB120');
y_limits = ylim();
axis tight;
ylim([0.8 * optimal_social_cost, y_limits(2)]);
axes = gca;
axes.XTickLabel = {'NE', 'Best possible'};
axes.Title.FontName = 'ubuntu';
axes.Title.String = 'Social cost';
axes.Title.FontSize = 20;
axes.XAxis.FontSize = 20;
axes.XAxis.FontName = 'ubuntu';
axes.YAxis.FontSize = 20;
axes.YAxis.FontName = 'ubuntu';