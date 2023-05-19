function [p_win, p_yield] = get_payments(b_win, b_yield)
% Payment function
%   Inputs:
%       b_win: winning bid
%       b_yield: losing bid
%   Outputs:
%       p_win: payment of winning player. Must be <= b_win, integer
%       p_yield: payment of yielding player. Must be <= b_yield, integer
%   Example: winning player pays its bid to yielding player (default case)
%       p_win = b_win;
%       p_yield = -b_win;
%   Instructions:
%       - For peer to peer payment set p_yield = -p_win.
%       - Some payment can be made to all the players in the game (instead
%           of just the yielding player).
%           This is achieved by letting p_yield != -p_win.
%           The yielding player will receive (-p_yield).
%           The surplus (p_win + p_yield) gets redistributed to all the
%           players in the game. (p_win + p_yield) must be >= 0.

%% Default: pay bid to peer - TO BE REPLACED WITH YOUR PAYMENT RULE
 p_win = b_win
 p_yield = 0;
 if(b_win > 7)
     p_win = floor(b_win / 2);
     
 end

end