function U = getU(obj, t)
if t <= obj.t1
    U = obj.U_max;
elseif t <= obj.t1 + obj.t2
    U = obj.U_max*((obj.t1 + obj.t2 - t)/obj.t2);
else
    U = 0;
end
end