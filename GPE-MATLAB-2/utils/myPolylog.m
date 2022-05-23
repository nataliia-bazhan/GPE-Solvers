function res = myPolylog(s, z)

N = size(z, 1);
r = linspace(-1, 1, N);
res = zeros('like', z);
for i = 1 : 100
    res = res + z.^i/i^s;
end

