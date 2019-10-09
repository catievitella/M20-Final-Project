function hE = create_hE(n,qnew)
    %creates energy hessian for any number of nodes
    
    %important constants
    l = 0.10; %total length of beam (m) 
    dl = l/(n-1); %length of each partition of beam (m) 
    r0 = 0.001; %cross-sectional radius of beam (m)
    E = 1*(10^9); %Young modulus (Pa)
    EA = E*pi*(r0)^2; %stretching stiffness
    EI = E*pi*((r0)^4)/4; %bending stiffness
    
    %hess
        hE = zeros(n*2);
        hEs = zeros(n*2);
        hEb = zeros(n*2);

        for k = 1:(n*2)
            if (mod(k,2) == 0) && (k >= 4)
                hEs(k-3:k, k-3:k) = hEs(k-3:k, k-3:k) + hessEs(qnew(k-3), qnew(k-2), qnew(k-1), qnew(k), dl, EA);
            end
            if (mod(k,2) == 0) && (k >= 6)
                hEb(k-5:k, k-5:k) = hEb(k-5:k, k-5:k) + hessEb(qnew(k-5), qnew(k-4), qnew(k-3), qnew(k-2), qnew(k-1), qnew(k), dl, EI);
            end
        end
        hE = hEs + hEb;
end