syms a b c d e w
assume(a ~= 0 & b ~= 0 & c ~= 0 & d ~= 0 & e ~= 0 & w ~= 0)

pol = a^2*(1+w^2)^2 + 2*a*c*(1+w^2)*w + c^2*w^2 ...
    - 4*(1-2*e)*(w^2 + (b-2*d)/(1-2*e)*w + 1) * (e*(1+w^2) + d*w);

pol_w4 = collect(pol, w);
simplify(subs(pol_w4, w, 1))

pol_d2 = pol_w4 / w^2;
