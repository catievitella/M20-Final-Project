function dE = create_dE(n, qnew)
    %creates energy gradient for any number of nodes
    %Catie Vitella
    %UID: 8050176562
    
    %important constants
    l = 0.10; %total length of beam (m) 
    dl = l/(n-1); %length of each partition of beam (m) 
    r0 = 0.001; %cross-sectional radius of beam (m)
    E = 1*(10^9); %Young modulus (Pa)
    EA = E*pi*(r0)^2; %stretching stiffness
    EI = E*pi*((r0)^4)/4; %bending stiffness
    
    %grad
    dE = zeros((n*2),1);
    dEs = zeros((n*2),1);
    dEb = zeros((n*2),1);
    
    for k = 1:(n*2)
            if (mod(k,2) == 0) && (k >= 4)
                dEs(k-3:k,1) = dEs(k-3:k,1) + gradEs(qnew(k-3), qnew(k-2), qnew(k-1), qnew(k), dl, EA);
            end
            if (mod(k,2) == 0) && (k >= 6)
                dEb(k-5:k,1) = dEb(k-5:k,1) + gradEb(qnew(k-5), qnew(k-4), qnew(k-3), qnew(k-2), qnew(k-1), qnew(k), dl, EI);
            end
        end
        dE = dEs + dEb;
end