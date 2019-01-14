function vv = RepeatArrayElements(v,rep)
% repeats the elements of v according to the values in rep

    v(rep==0) = [];
    rep(rep==0) = [];
    index = zeros(1,sum(rep));
    index(cumsum([1 rep(1:end-1)])) = 1;
    index = cumsum(index);
    vv = v(index);
end