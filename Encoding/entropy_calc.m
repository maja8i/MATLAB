%-------------------------------------------------------------------------%

%Function for computing theoretical (practical) entropy of a given set of 
%input values, in bits.

%---- INPUTS ----

%input: Vector of input values that we wish to compute entropy for.

%---- OUTPUTS ----

%entropy_out: Computed entropy, in bits. This is the average minimum number
%             of bits needed PER SYMBOL.

%-------------------------------------------------------------------------%

function [entropy_out] = entropy_calc(input)

%Initialize "entropy_out" to 0
entropy_out = 0;

%Extract the unique symbols in the input vector
symbols = unique(input);
%Compute entropy
for i = 1:length(symbols)
    p = (numel(find(input == symbols(i))))/numel(input);    %'Probability' of getting the current value
    entropy_out = entropy_out + p*(-log2(p));
end
