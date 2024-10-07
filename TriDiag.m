function A = TriDiag(n)
    % 创建n*n三对角矩阵
    mainDiag = 4 * ones(n, 1);  % 主对角线上的元素
    offDiag = -1 * ones(n - 1, 1);  % 副对角线上的元素
    T = diag(mainDiag) + diag(offDiag, 1) + diag(offDiag, -1); 

    % 创建 n x n 单位矩阵
    I = -eye(n); 

    % 初始化整个?n^2 x n^2 矩阵
    A = zeros(n^2);

    % 填充矩阵
    for i = 1:n
        A((i-1)*n+1:i*n, (i-1)*n+1:i*n) = T;
        if i < n
            A(i*n+1:(i+1)*n, (i-1)*n+1:i*n) = I;
            A((i-1)*n+1:i*n, i*n+1:(i+1)*n) = I;
        end
    end
end
