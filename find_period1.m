function [T] = find_period1(seq,TT)
i = 1;
while seq(i) < seq(i+1) 
    i = i+1;
    if i>length(seq)-1
        break
    end
end
i
T = 2*TT(i);

    
end