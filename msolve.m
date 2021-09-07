function [F] = msolve(x,TR,KBN,KSN,KCN,KMN,TNa,TBa,TSr,TCa,TMg)

F = [x(1)+x(2)+x(3)+x(4)+x(5)-TR;
x(2)*x(6)^2-KBN*x(1)^2*x(7);
x(3)*x(6)^2-KSN*x(1)^2*x(8);
x(4)*x(6)^2-KCN*x(1)^2*x(9);
x(5)*x(6)^2-KMN*x(1)^2*x(10);
x(1)+x(6)-TNa;
x(2)+x(7)-TBa;
x(3)+x(8)-TSr;
x(4)+x(9)-TCa;
x(5)+x(10)-TMg];

end