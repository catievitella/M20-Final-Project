function R = createR(dl,nodes)
%creates R vector for any number of nodes
    mid = (nodes/2)+0.5;
    R = zeros(nodes,1);
    R(1:mid-1) = dl/10;
    R(mid) = 0.025;
    R(mid+1:nodes) = dl/10;
end