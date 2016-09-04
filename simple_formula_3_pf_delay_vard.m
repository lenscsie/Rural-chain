% calculate probability of delivery for a simple chain using the formula
% Prob = 0.5*(1 - erf((10*log10((2^R - 1)*sigma_awgn) - alpha + beta.*log10([1:Distance]))/(sqrt(2)*sigma))); % Ps = 1 - P_outage
I = find(Prob >= Prob_th);

for i = N_max:-1:1%length(d)
    if isempty(I) ~= 1
        d = I(length(I))/i; % Internodal distance
        n = ceil(Distance/d); % No of nodes
         
        
        % T = range where a transmission could be consider received, depending on
        % Prob_th
    
        P_dest_k = zeros(1,n*T_R); % Reception prob at the destination for each step 
    
        T = i;%I(length(I)); % T is in terms of number of nodes - range of transmission
        
        if flag_mimo == 1
            P = 0.5*(1 - erf((10*log10((2^(R/Nmin) - 1)*sigma_awgn*Nmin/Nrecv) - alpha + beta.*log10(d*[1:T]))/(sqrt(2)*sigma)));%Prob(I); % Probability of reception at nodes +1,..., +T
        else
            P = 0.5*(1 - erf((10*log10((2^R - 1)*sigma_awgn) - alpha + beta.*log10(d*[1:T]))/(sqrt(2)*sigma)));%Prob(I); % Probability of reception at nodes +1,..., +T
        end
        P = P*(1-pf);% Random node failures

        clear X0 X A P_d
        N = n;%ceil(Distance/d(i)); % No. of nodes to destinatio n
        K = ceil(N/T); % No. of iterations we need to get to the destination
        step_time = packet_size/bit_rate + min(T,N - [1:N])*ACK_size/bit_rate;
        
        X0(1) = prod(1 - P); % X0 - probability of being potential forwarder for nodes +1,...,+T from the transmiter, 
        for j = 1:T-1
            X0(j+1) = P(j) * prod(1 - P(j+1:T)); % j received it, and j+1:T did not receive it
        end
        X0(T+1) = P(T);
        X = [1; zeros(N-1,1)];%[X0'; zeros(N-length(X0),1)]; % X is a column vector  % X - probability of being potential forwarder,
        X_save = [X' zeros(1,N-length(X))]; % Full vector of 1:N nodes - prob of being a forwarder
        
        % Construction of A - reception matrix
        % A - N by N matrix
     
        clear A P_d
        %A(1,:) = [X0 zeros(1,N-T-1)];
        %A(N,:) = [zeros(1,N-1) (1-P(1))]; % We should consider that pkts are not reaching the destination
        P_d = zeros(1,N);
        %P_d(N) = P(1); % Probability of reaching the destination node N+1
        for j = N-1:-1:N-T
            tmp = prod(1 - P(1:N-j)); % We should consider that pkts are not reaching the destination
            for l = 1:N-j-1
                tmp = [tmp prod(1 - P(l+1:N-j))];
            end
            %tmp = [tmp 1];    
            A(j+1,:) = [zeros(1,j) [1 P(1:N-j-1)].*tmp];
            if N-j <= T
                P_d(j+1) = P(N-j); % Reception probabilities at destination
            end
        end
        for j = N-T-1:-1:0
            A(j+1,:) = [zeros(1,j) X0 zeros(1,N-T-1-j)]; %[X(T - j + 1:T) zeros(1,T-j)];
        end
        
        A = A';
        multiplier = diag(A); % we will need to limit re-transmissions
       
        %%%%%%%%%%
        P_dest_k(1,1) = P_d*X;
        avg_delay_k(1,1) = (step_time.*P_d)*X;
        % Rounds of opportunistic routing
        for k = 1:n*T_R-1 % -1 to take of the transmission to destination which is considered by P_d
           
            X = A*X;
            
            
            % Taking into account the re-transmission limit T_R
            if k == T_R  
                X = X - X_save(k-T_R+1,1:length(X))'.*multiplier.^T_R;
            elseif k > T_R 
                X = X - ((A - diag(multiplier))*X_save(k-T_R,1:length(X))').*multiplier.^T_R;% 
            end
            % Taking care of the float precision error 
            I_tmp = find(X < 1e-10); X(I_tmp) = 0;
             X_save = [X_save; X' zeros(1,N-length(X))];
             P_dest_k(1,k+1) = P_d*X;
             avg_delay_k(1,k+1) = k*(step_time.*P_d)*X;
        end
        
%         P_n(i) = P_dest_k(1);
%         for j = 2:length(P_dest_k)
%             P_n(i) = P_n(i) + P_dest_k(j) * prod(1-P_dest_k(1:j-1));
%         end
         P_n(i,height_index,Power_index) = sum(P_dest_k);
         %step_time = packet_size/bit_rate + T*ACK_size/bit_rate;
         avg_delay(i,height_index,Power_index) = sum(avg_delay_k)/sum(P_dest_k);%P_dest_k.*[1:length(P_dest_k)])*step_time/sum(P_dest_k);
         n_save(i,height_index,Power_index) = n+1;
         d_save(i,height_index,Power_index) = d;
    else
        disp('Nodes are unreachable')
        break;
    end
        
end
