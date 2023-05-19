clear all
close all
clc

load('DSMgame_data.mat');

A = [];
for i = 1:100
    for j = 1:100
        if i == j
            B = 2*p1(i) * eye(24);
        else
            B = p1(i) * eye(24);
        end
        
        if i == 1 && j ==1
            C = B;
        else
            C = [A; B];
        end
    end
    
    if i == 1
        A = C;
    else
        A = [A C];
    end
end
