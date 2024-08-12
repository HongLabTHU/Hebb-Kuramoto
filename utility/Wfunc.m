function W = Wfunc(v, u)
a = abs(Aff_func(v, u));
b = abs(Aff_func(v, conj(u)));
W = a^2 - b^2;
end



