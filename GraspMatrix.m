function G = GraspMatrix(EP)
    G = [[eye(3),   zeros(3,3)];
        [Hat(EP),   eye(3)]];
end