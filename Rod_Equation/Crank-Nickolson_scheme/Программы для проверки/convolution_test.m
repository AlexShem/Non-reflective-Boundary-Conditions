clear;
p1 = [.1, .5, .9];
n = 30;

p_true = conv(p1, p1);
p_cut = p_true;
p_cut = p_cut(end-length(p1)+1:end);
for i = 3 : n
    p_true = conv(p_true, p1);
    p_cut = conv(p_cut, p1);
    p_cut = p_cut(end-length(p1)+1: end);
end
