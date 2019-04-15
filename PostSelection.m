function PostSelectOp = PostSelection(st, E)
    Z = {[1;0], [0;1]};
    U = PostSelectionU(E);
    PI = kron(eye(2), Z{1}*Z{1}');
    UPU = U'*PI*U;
    inner = UPU * kron(st, Z{1}*Z{1}') * UPU;
    PostSelectOp = (1/trace(E*st)) * PartialTrace(inner);
end

function PostUOp = PostSelectionU(MO)
    Z = {[1;0], [0;1]};
    PostUOp = zeros(4);
    for idx = [1,2]
        PostUOp = PostUOp + kron(SqrtMat(MO),eye(2)) * kron(Z{idx},Z{1}) * kron(Z{idx},Z{1})';
        PostUOp = PostUOp + kron(SqrtMat(eye(2) - MO), eye(2)) * kron(Z{idx},Z{2}) * kron(Z{idx},Z{1})';
        PostUOp = PostUOp + kron(Z{idx}, Z{2}) * kron(Z{idx}, Z{2})';
    end
end

function SqrtOp = SqrtMat(Operator)
    [V,~,~] = eig(Operator);
    es = eig(Operator);
    SqrtOp = [0,0;0,0];
    for idx = (1:numel(es))
        vec = V(:, idx);
        mat = vec * vec';
        SqrtOp = SqrtOp + sqrt(es(idx))*mat;
    end
end

