function p = polysum(u, v)
    n = abs(length(u) - length(v));
    if length(u) > length(v)
        p = [u(1:n), u(n+1:end)+v];
    else
        p = [v(1:n), v(n+1:end) + u];
    end
end