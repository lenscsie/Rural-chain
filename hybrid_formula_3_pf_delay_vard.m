% calculate probability of delivery for a hybrid chain using the formula
% Prob = 0.5*(1 - erf((10*log10((2^R - 1)*sigma_awgn) - alpha + beta.*log10([1:Distance]))/(sqrt(2)*sigma))); % Ps = 1 - P_outage
I = find(Prob >= Prob_th);

for i = N_max:-1:1%1:length(d)
    if isempty(I) ~= 1
        d = I(length(I))/i; % Internodal distance
        n = ceil(Distance/d); % No of nodes
        
        clear X0 X A P_d
            if M_hybrid > 0
                N = n + floor((n-1)/M_hybrid);% No. of nodes to destination, destination doesn't have a redundant peer
                step_time = packet_size/bit_rate + ...
                (min(T,N-[1:N]) + sum(rem(kron([1:N]',ones(1,T))+kron([1:T],ones(N,1)),M_hybrid)== 0,2)')*ACK_size/bit_rate;
        
            else
                N = n;
                step_time = packet_size/bit_rate + min(T,N - [1:N])*ACK_size/bit_rate;
            end
        
        % T = range where a transmission could be consider received, depending on
        % Prob_th
        P_dest_k = zeros(1,N*T_R); % Reception prob at the destination for each step
    
        T = i;%I(length(I)); % T is in terms of number of nodes - range of transmission
        if flag_mimo == 1
            P = 0.5*(1 - erf((10*log10((2^(R/Nmin) - 1)*sigma_awgn*Nmin/Nrecv) - alpha + beta.*log10(d*[1:T]))/(sqrt(2)*sigma)));%Prob(I); % Probability of reception at nodes +1,..., +T
        else
            P = 0.5*(1 - erf((10*log10((2^R - 1)*sigma_awgn) - alpha + beta.*log10(d*[1:T]))/(sqrt(2)*sigma)));%Prob(I); % Probability of reception at nodes +1,..., +T
        end
        
        % For redundant peers at mth position
        Prob_redundant = 0.5*(1 - erf((10*log10((2^R - 1)*sigma_awgn) - alpha + beta.*log10(((d.*[1:T]).^2 + d_peer^2).^(0.5)))/(sqrt(2)*sigma))); % Ps = 1 - P_outage
        Prob_peer = 0.5*(1 - erf((10*log10((2^R - 1)*sigma_awgn) - alpha + beta.*log10(d_peer))/(sqrt(2)*sigma)));

        P_r = Prob_redundant([1:T]); % For redundant peers
        P = P*(1-pf);% Random node failures
        P_r = P_r*(1-pf);
        Prob_peer = Prob_peer*(1-pf);

        
        %K = ceil(N/T); % No. of iterations we need to get to the destination

        X = [1; zeros(N-1,1)]; % X is a column vector  % X - probability of being potential forwarder,
        X_save = [X' zeros(1,N-length(X))]; % Full vector of 1:N nodes - prob of being a forwarder

        % Construction of A - reception matrix
        % A - N by N matrix

%         A(N,:) = [zeros(1,N-1) (1-P(1))]; % We should consider that pkts are not reaching the destination
%         if rem(n-1,M_hybrid) == 0
%             A(N-1,:) = [zeros(1,N-2) (1-P_r(1))*(1-Prob_peer) Prob_peer*(1-P_r(1))];
%         end
% 
         P_d = zeros(1,N);
%         P_d(N) = P(1); % Probability of reaching the destination node N+1
%         if rem(n-1,M_hybrid) == 0
%             P_d(N-1) = P_r(1);
%         end
        
        for j = 0:n-1%-2
            tlim = min([n-j, T]);% We should consider that pkts are not reaching the destination
            tlim2 = min([n-j-1, T]); % Not a good way to handle the issue, destination does not have a redundant peer
            clear X0
            tmp = zeros(1,tlim2); % Do we have redundant peer within the range or not?
            for l = 1:length(tmp)
                if rem(j+l,M_hybrid) == 0
                    tmp(l) = P_r(l);
                end
            end
            
            X0 = prod(1 - P(1:tlim))*prod(1 - tmp); 
            for l = 1:tlim-1
                if rem(j+l,M_hybrid) == 0
                    X0 = [X0 tmp(l)* prod(1 - P(l:tlim))*prod(1 - tmp(l+1:tlim2))];
                end
                X0 = [X0 P(l)*prod(1 - P(l+1:tlim))*prod(1 - tmp(l+1:tlim2))];
            end
            if j < n - T
                if rem(j+T,M_hybrid) == 0
                    X0 = [X0 tmp(T)*(1-P(T))];
                end
                X0 = [X0 P(T)];
            end
            if M_hybrid > 0
                index = j+floor(j/M_hybrid);
            else
                index = j;
            end
            A(index+1,:) = [zeros(1,index) X0 zeros(1,N-length(X0)-index)];
            if n-j <= T
               
                P_d(index+1) = P(n-j); % Reception probabilities at destination
            end
            
            if (rem(j,M_hybrid) == 0 && j ~= 0) % Redundant peer will transmit
                clear X0
                tmp = zeros(1,tlim2);
                for l = 1:tlim2
                    if rem(j+l,M_hybrid) == 0
                        tmp(l) = P(l);
                    end
                end

                X0(1,1) = 1* prod(1 - P_r(1:tlim))*prod(1 - tmp)*(1-Prob_peer); % X0 - probability of being potential forwarder for nodes +1,...,+T from the transmiter,
                X0(1,2) = Prob_peer*prod(1 - P_r(1:tlim))*prod(1-tmp);
                for l = 1:tlim-1
                    if rem(j+l,M_hybrid) == 0
                        X0 = [X0 tmp(l)* prod(1 - P_r(l:tlim))*prod(1 - tmp(l+1:tlim2))];
                    end
                    X0 = [X0 P_r(l) * prod(1 - P_r(l+1:tlim))* prod(1 - tmp(l+1:tlim2))]; % j received it, and j+1:T did not receive it

                end
                if j < n - T
                    if rem(j+T,M_hybrid) == 0
                        X0 = [X0 tmp(T)*(1-P_r(T))];
                    end
                    X0 = [X0 P_r(T)];
                end
                if M_hybrid > 0
                    index = j+floor(abs(j-1)/M_hybrid);
                else
                    index = j;
                end
                A(index+1,:) = [zeros(1,index) X0 zeros(1,N-length(X0)-index)];
                if n - j <= T
                    P_d(index + 1) = P_r(n - j); % Reception at destination
                end
            end

            
        end

        A = A';
        multiplier = diag(A); % we will need to limit re-transmissions

        %%%%%%%%%%
        P_dest_k(1,1) = P_d*X;
        avg_delay_k(1,1) = (step_time.*P_d)*X;
        % Rounds of opportunistic routing
        for k = 1:N*T_R-1 % -1 to take of the transmission to destination which is considered by P_d

            X = A*X;


            % Taking into account the re-transmission limit T_R
            if k == T_R
                X = X - X_save(k-T_R+1,1:length(X))'.*multiplier.^T_R;
            elseif k > T_R
                X = X - ((A - diag(multiplier))*X_save(k-T_R,1:length(X))').*multiplier.^T_R;%
            end
            % Taking care of the float precision error
            %I_tmp = find(X < 1e-10); X(I_tmp) = 0;
            X_save = [X_save; X' zeros(1,N-length(X))];
            P_dest_k(1,k+1) = P_d*X;
            avg_delay_k(1,k+1) = k*(step_time.*P_d)*X;
        end

        %         P_n(i) = P_dest_k(1);
        %         for j = 2:length(P_dest_k)
        %             P_n(i) = P_n(i) + P_dest_k(j) * prod(1-P_dest_k(1:j-1));
        %         end
        P_n_hybrid(i,height_index,Power_index) = sum(P_dest_k);
        %step_time = packet_size/bit_rate + (T + ceil(T/M_hybrid))*ACK_size/bit_rate;
        avg_delay_hybrid(i,height_index,Power_index) = sum(avg_delay_k)/sum(P_dest_k);%PP_dest_k.*[1:length(P_dest_k)])*step_time;%/sum(P_dest_k);
    else
        disp('Nodes are unreachable')
        break;
    end

end
