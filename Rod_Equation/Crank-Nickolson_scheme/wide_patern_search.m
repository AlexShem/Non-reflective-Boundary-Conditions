layers = 4;
start = 1;
time = 6;

for i = start : time+1
    for j = start : time+1
        for k = start : time+1
            for s = start : time+1
                deg = [i-1, j-1, k-1, s-1];
                coeffs_contour_plot(deg, deg);
                disp(['Is Done: [' num2str(deg(1)), ', ' num2str(deg(2)) ', '...
                    num2str(deg(3)) ', ' num2str(deg(4)) ']']);
            end 
        end 
    end 
end