function f = pm_constraint_vector(X, R)
%pm_vector_constraints Get PM constraints on R, X as a vector function
n = size(X, 1);
d = size(R, 1);
f = [];

for idx=1:n
    f = [f; sum(X(:,idx))-1; sum(X(idx,:))-1];
    for jdx=1:n
        f = [f; X(idx,jdx)*(X(idx,jdx) - 1)];
    end
end

for idx=1:d
    f = [f; R(:,idx).'*R(:,idx) - 1];
    for jdx=1:idx
        f = [f; R(:,idx).'*R(:,jdx)];
    end
end

end

