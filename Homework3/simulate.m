clear all
close all
clc


%Initialization
values = zeros(21,1);
for i  = 0:20
    values(i + 1) = (23 / 5 * (41 - 2*i))/42;
end



optimal_high_bid = zeros(21,1);
optimal_low_bid  = zeros(21,1);

%Backward induction
steps = 19;

for s = steps:-1:1
    new_values = zeros(21,1);
    for k_i = 0:20
        values_high = zeros(k_i + 1,1);
        for b_i = 0:k_i
            values_high(b_i + 1) = (41-2*b_i)*10/42 + dot(TransitionProbabilities(k_i,b_i),values);
        end
        values_low = zeros(k_i + 1,1);
        for b_i = 0:k_i
            values_low(b_i + 1) = (41-2*b_i)/42 + dot(TransitionProbabilities(k_i,b_i),values);
        end
       
        new_values(k_i + 1) = 0.4 * min(values_high) + 0.6 * min(values_low);
         
        %Printing optimal strategy at the beginning
        if s == 1
            [optimal_high_value,optimal_high_bid_idx] = min(values_high);
            [optimal_low_value,optimal_low_bid_idx] = min(values_low);
            optimal_high_bid(k_i + 1) = optimal_high_bid_idx - 1;
            optimal_low_bid(k_i + 1) = optimal_low_bid_idx - 1;
        end
        
    end
    values = new_values
end

optimal_high_bid
optimal_low_bid

%var = TransitionProbabilities(16,11)

%Initial karma level of 100 cars
karma = 10*ones(100,1);

%Defining number of iterations
time_steps = 1000;

%Initializing an array to kip track of the bids
bids_tracker = [];

%Initializing a matrix keeping count of each player cost in the last 100
%iterations, matrix is 100 x 100 ( Players X scenarios)
players_costs = zeros(100,100);

%Initializing a counter keeping track of the cumulative cost
social_cost = 0;

%External loop over the time steps
for timestep = 1:time_steps
    
    %Candidates is a vector containing the indexes of the 100 cars. It is used
    %to sample the couples of players.
    candidates = linspace(1,100)';

    for couple_index = 1:50

        %Sampling a random couple
        couple = randsample (candidates,2);
        index_player1 = couple(1);
        index_player2 = couple(2);
        
        %Random process to determine the urgency of the 2 players
        urgency1 = randsample(10,1);
        urgency2 = randsample(10,1);

        %Retrieving karma values of the two players
        k1 = karma(index_player1);
        k2 = karma(index_player2);
        
        %Choosing Player 1 optimal action according to his karma value
        if urgency1 <= 6
            b1 = optimal_low_bid(k1 + 1);
        else
            b1 = optimal_high_bid(k1 + 1);
        end
        
        %Choosing Player 2 optimal action according to his karma value
        if urgency2 <= 6
            b2 = optimal_low_bid(k2 + 1);
        else
            b2 = optimal_high_bid(k2 + 1);
        end
        
        %Keeping track of the bids after the 900 iteration
        if timestep > time_steps - 100
            bids_tracker = [bids_tracker; b1; b2];
        end
        %Retrieving P1 and P2 new karma values
        [newk1, newk2, winner] = ComputenewKarma(k1,b1,k2,b2);

        %Computing cost for this interaction
        
         if timestep > time_steps - 100
             [cost1,cost2] = ComputeCost(urgency1, urgency2, winner);
             players_costs(timestep - time_steps + 100,index_player1) = cost1;
             players_costs(timestep - time_steps + 100, index_player2) = cost2;
             normalized_cost = (cost1 + cost2)/ 100;
             social_cost  = social_cost + normalized_cost;
        end
        
        %Updating karma vector
        karma(index_player1) = newk1;
        karma(index_player2) = newk2;

        candidates = setdiff(candidates,couple);
    end
    karma = capAndRedistribute(karma,20);
    
end
SC = social_cost / 100
sum_player_cost = sum(players_costs,1);
max_cost = max(sum_player_cost)
min_cost = min(sum_player_cost)
histogram(bids_tracker);
x = [0:1:20]
figure
scatter(x,optimal_high_bid)
hold on
scatter(x,optimal_low_bid)


function [cost1,cost2] = ComputeCost(urgency1, urgency2, winner)
    if winner == 1
        if urgency2 > 6
            cost1 = 0;
            cost2 = 10;

        else
            cost1 = 0;
            cost2 = 1;
        end
    elseif winner == 2
        if urgency1 > 6
            cost1 = 10;
            cost2 = 0;
        else
            cost1 = 1;
            cost2 = 0;
        end
    end
end
    


    function [newk1,newk2, winner] = ComputenewKarma(k1,b1,k2,b2)
%Given the karma valeus and the bids of the two players, this function
%return determines the winner of the auction and computes the new karma
%without capping it.

    if(b1 > b2)
        newk1 = k1 - b1;
        newk2 = k2 + b1;
        winner = 1;
    elseif (b1 == b2)
            %Tossing a coin between 1(heads) and 2(tails) to break the tie. If
            % heads Player 1 wins, if tails Player 2 wins.
            coin_toss = randsample (2,1);
            if coin_toss == 1
                newk1 = k1 - b1;
                newk2 = k2 + b1;
                winner = 1;
            else
                newk1 = k1 + b2;
                newk2 = k2 - b2;
                winner = 2;
            end
    else
        newk1 = k1 + b2;
        newk2 = k2 - b2;
        winner = 2;
    end
end



function probs = TransitionProbabilities (k_i,b_i)
        if(b_i > k_i)
            fprintf("RAISE VALUE ERROR,Cannot bet more than you have")
        end

        % Given k_i and b_i returns a vector of transition probabilities
        % P(k_i / K_i, b_i)
        probs = zeros(21,1);
        if k_i + b_i >= 20
            if(b_i == 0)
                probs(21) = 1;
            else
            probs(k_i - b_i + 1) = (1 + 2 * b_i)/42;
            probs(21) = (41 - 2 * b_i) / 42;
            end
        
        else
            if(b_i == 0)
                probs(k_i + 1) = (1 + 2 * b_i)/42 + 1 / 42;
            else

                probs(k_i - b_i + 1) = (1 + 2 * b_i)/42;
                probs(k_i + b_i + 1) = 1 / 42;
            end
            index = k_i + b_i + 2;
            while index < 21
                probs(index) = 1 / 21;
                index = index + 1;
            end
            probs(21) = (k_i + 1) / 21;

  
        
        end
        

end



function newk = capAndRedistribute(k, maxvalue)

% capAndRedistribute receives a vector of length N with the karma of % all players, caps the karma at maxvalue, and redistributes the
% excess to the other players. The function returns the resulting
% vector of karma counters.

    newk = min(k,maxvalue);
    excess = sum(k) - sum(newk);
    while excess > 0
        recipients = find(newk < maxvalue);
        if excess >= length(recipients)
            newk(recipients) = newk(recipients) + 1; 
            excess = excess - length(recipients);
        else
            recipients = recipients(randperm(length(recipients))); 
            newk(recipients(1:excess)) = newk(recipients(1:excess)) + 1; 
            excess = 0;
        end
    end
end