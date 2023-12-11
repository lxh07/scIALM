function [I J col omega] = allsample(m, n, p, mask) 

[row,col]=find(mask);
I = row;
J = col;

omega = sub2ind([m,n],I,J);
col = [0; find(diff(J)); p];   %diff():求差分或求导。find():寻找J中非零元素，并返回其索引值（按列排序）。