for i = 1:10 
    fprintf("%d", i);
end

function y=fun1(x)
  y=x;
end

[x1] = fun1(3);
[x2] = fun1(3);
[x3] = fun1(3);