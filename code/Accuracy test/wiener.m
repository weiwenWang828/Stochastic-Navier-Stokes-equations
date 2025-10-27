function [wiener1,wiener2] = wiener(xx,yy,k0,N,params)

xi = randn(params.J, params.J);
w1 = zeros(N, N);
w2 = w1;


for i = 1 : params.J
    for j = 1 : params.J
        w1 = w1 + xi(i, j) .* sin(i * pi * (xx + 1) / 2) .* sin(j * pi * (yy + 1) / 2) / (i + j) .^ sqrt(2 + params.epsilon); 
        w2 = w2 + xi(i, j) .* sin(i * pi * (xx + 1) / 2) .* sin(j * pi * (yy + 1) / 2) / (i + j) .^ sqrt(2 + params.epsilon); 

    end
end

wiener1 = sqrt(k0) * w1;
wiener2 = sqrt(k0) * w2;

