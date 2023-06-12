
function out = make_nan(in)


 I = find(in);    out = nan(size(in));  out((I)) = 1;


