function plot_line(x,v,mylimits_x,mycolor)

w = null(v);
P = mylimits_x;
X = x(1) + w(1)*P;
Y = x(2) + w(2)*P;
plot(X,Y,mycolor);

end

