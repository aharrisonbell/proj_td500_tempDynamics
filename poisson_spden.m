% y=g1dfun(x,Xmax,A,Sigma)
%    One dimensional poisson function. 
%    Based on Thompson and Schall
%    Xmax  = Position of the Peak
%    A     = Peak maximum

function y = poisson_spden(time, Spike);
y = zeros(1,length(time));
for t = Spike:length(time),
    B = Spike - time(t) - 1;
    A = B/20;
    y(t) = (1 - exp(B)) * exp(A) / 19.05;
end
return

