function y = UVtOmega(M,I,J,col);

y = zeros(length(I), 1);
for k = 1:length(col)-1
    j = J(col(k)+1);
    Xj = M(:,j);
    idx = [col(k)+1:col(k+1)];
    y(idx) = Xj(I(idx));
end
