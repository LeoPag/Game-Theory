function tax = get_karma_tax(karma)
% Payment function
%   Inputs:
%       karma: karma level of player
%   Outputs:
%       tax: tax to deduct from player and redistribute to all the players
%            in the game. Must be <= k
%   Example: linear tax with rate 0.1
%       tax = 0.1 * karma;
%   NOTE: A fractional tax will be interpolated probabilstically.
%         Example: tax = 2.3 means that 2 will be deducted with
%         probability 0.7 and 3 with probability 0.3

%% Default: no taxation - TO BE REPLACED WITH YOUR TAX RULE
tax = 0.1 * karma;


