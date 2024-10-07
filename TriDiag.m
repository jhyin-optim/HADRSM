function A = TriDiag(n)
    % ����n*n���ԽǾ���
    mainDiag = 4 * ones(n, 1);  % ���Խ����ϵ�Ԫ��
    offDiag = -1 * ones(n - 1, 1);  % ���Խ����ϵ�Ԫ��
    T = diag(mainDiag) + diag(offDiag, 1) + diag(offDiag, -1); 

    % ���� n x n ��λ����
    I = -eye(n); 

    % ��ʼ������?n^2 x n^2 ����
    A = zeros(n^2);

    % ������
    for i = 1:n
        A((i-1)*n+1:i*n, (i-1)*n+1:i*n) = T;
        if i < n
            A(i*n+1:(i+1)*n, (i-1)*n+1:i*n) = I;
            A((i-1)*n+1:i*n, i*n+1:(i+1)*n) = I;
        end
    end
end
