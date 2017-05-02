function c = inner(a,b)
% Input: The two vectors a and b 
% Output: The value of the inner product of a and b.

c = nansum(a.*b);
end