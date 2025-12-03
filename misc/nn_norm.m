function val = nn_norm(t, x1, x2)
    % -- val = NN_NORM(t, x1, x2)
    %   Calculates the nearest neighbor norm of x1 and x2
    N = min([length(t), size(x1, 1), size(x2, 1)]);

    %calculate NN metric
    val = 0;
    for i = 1:N
        % calculate distance from one point to another and taking the
        % distance to the nearest one as a part of norm
        val = val + min(vecnorm(x1(i, :) - x2, 2, 2));
    end
    val = val / N; % mean by N
end

