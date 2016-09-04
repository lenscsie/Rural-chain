% Calculates empirical delivery probability for simple chain
% Effective capacity in saturation condition

% Prob = 0.5*(1 - erf((10*log10((2^R - 1)*sigma_awgn) - alpha + beta.*log10([1:Distance]))/(sqrt(2)*sigma))); % Ps = 1 - P_outage
I = find(Prob >= Prob_th);



for i = N_max:-1:1%1:length(d)
    %ndelivered = 0;

    if isempty(I) ~= 1
        d = I(length(I))/i; % Internodal distance
        n = ceil(Distance/d); % No of nodes
        T = i;%I(length(I)); % T is in terms of number of nodes - range of transmission
        
         if flag_mimo == 1
            P = 0.5*(1 - erf((10*log10((2^(R/Nmin) - 1)*sigma_awgn*Nmin/Nrecv) - alpha + beta.*log10(d*[1:T]))/(sqrt(2)*sigma)));%Prob(I); % Probability of reception at nodes +1,..., +T
        else
            P = 0.5*(1 - erf((10*log10((2^R - 1)*sigma_awgn) - alpha + beta.*log10(d*[1:T]))/(sqrt(2)*sigma)));%Prob(I); % Probability of reception at nodes +1,..., +T
         end
        
        % General case, all possible signals
        capacity_tmp = zeros(1,Ntrials); % reached packets/time taken
        for j = 1:Ntrials
            delivered = zeros(1,n); % status of nodes having packet or not, 0 if forwarded
            delivered(1,1) = 1;
            forwarder = 1;
            t_count = zeros(1,n); % No. of trans from each node
            nreached = 0; %no. of packets reaching the destination
            ndropped = 0; %no. of packets dropped
            steps = 0;
            while(nreached+ndropped < Npackets)
                % forwader list for this round
                forwarder = find(delivered);
                % re-transmission consideration
                for k = 1:length(forwarder)
                    if t_count(forwarder(k)) >= T_R % drop the packet
                        ndropped = ndropped + 1;
                        if forwarder(k) ~= 1
                            delivered(1,forwarder(k)) = 0;
                        end
                        t_count(forwarder(k)) = 0;
                    end
                end
                forwarder = find(delivered); % updated list
                % check who can transmit without interference, prefrence
                % will be given to farther away nodes
                tmp = forwarder(length(forwarder));
                for k = length(forwarder)-1:-1:1
                    if forwarder(k+1)-forwarder(k) > 2*T % interferring range
                        tmp = [tmp forwarder(k)];
                    end
                end
                forwarder = tmp;
                clear tmp

                for k = 1:length(forwarder)
                    M = min([T   n+1-forwarder(k)]);
                    tmp = ((rand(1,M) <= P(1:M)).*(rand(1,M) <= (1-pf)));
                    t_count(forwarder(k)) = t_count(forwarder(k)) + 1;
                    I_del = find(tmp);
                    if isempty(I_del) ~= 1
                        if forwarder(k) ~= 1 % source always has a packet to send
                            delivered(1,forwarder(k)) = 0; % packet has been forwarded
                        end
                        t_count(forwarder(k)) = 0; % reset no. of transmissions
                        if forwarder(k) + I_del(length(I_del)) == n + 1
                            nreached = nreached + 1;
                        else
                            delivered(1,forwarder(k) + I_del(length(I_del))) = 1;
                        end
                    end
                end
                if nreached >= 1 % just to cut of the time taken for first packet.
                    steps = steps+1;
                end

            end
            ndelivered(j) = nreached;%ndelivered = ndelivered + nreached;
            if nreached == 0
                capacity_tmp(j) = 0;
            else
                capacity_tmp(j) = nreached/steps;
            end
        end

        % To calculate capacity in terms of bit/s
        % time taken in each step
        step_time = packet_size/bit_rate + T*ACK_size/bit_rate;
        capacity(i,height_index,Power_index) = mean(capacity_tmp*packet_size/step_time*1e-6); % in Mbits/s
        capacity_std(i,height_index,Power_index) = std(capacity_tmp*packet_size/step_time*1e-6);
        P_emp(i,height_index,Power_index) = mean(ndelivered/Npackets);%ndelivered/Ntrials/Npackets;
        P_emp_std(i,height_index,Power_index) = std(ndelivered/Npackets);
        
    end
end

