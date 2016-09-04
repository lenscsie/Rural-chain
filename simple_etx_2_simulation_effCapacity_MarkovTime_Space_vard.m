% Traditional ETX based routing
% We are not calculating ETX, rather using the probabilities of links,
% it is logically the same
% Prob = 0.5*(1 - erf((10*log10((2^R - 1)*sigma_awgn) - alpha + beta.*log10([1:Distance]))/(sqrt(2)*sigma))); % Ps = 1 - P_outage
I_p = find(Prob >= Prob_th);

for i = N_max:-1:1%1:length(d);
    %ndelivered = 0;
    d = I_p(length(I_p))/i; % Internodal distance
    n = ceil(Distance/d); % No of nodes
        
    if flag_mimo == 1
        Prob_tmp = 0.5*(1 - erf((10*log10((2^(R/Nmin) - 1)*sigma_awgn*Nmin/Nrecv) - alpha + beta.*log10(d*[1:n]))/(sqrt(2)*sigma)));%Prob(I); % Probability of reception at nodes +1,..., +n
    else
        Prob_tmp = 0.5*(1 - erf((10*log10((2^R - 1)*sigma_awgn) - alpha + beta.*log10(d*[1:n]))/(sqrt(2)*sigma)));%Prob(I); % Probability of reception at nodes +1,..., +n
    end
    
    I = find(Prob_tmp > 1e-3);
    if isempty(I) ~= 1
        T = I(length(I));
        Prob_tmp = Prob_tmp(I);

        % To determine the best ETX path
        hop1 = [1 1];

        for j = 2:T
            if Prob_tmp(j) > prod(Prob_tmp(hop1))
                hop1 = j;

            end
            if j < T
                hop1 = [hop1 1];
            end
            % No need to fix the tail here, as we are only interested in the
            % max. hop with highest probability i.e. max(hop1);
        end

        hop = max(hop1);


        %%%%% Simulation start %%%%%%%
        capacity_tmp = zeros(1,Ntrials); % reached packets/time taken
        for j = 1:Ntrials
            delivered = zeros(1,n); % status of nodes having packet or not
            delivered(1,1) = 1;
            state = 2*ones(1,n+1); % node is failed or not 1 - failed, 2- not failed
          
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
                            delivered(forwarder(k)) = delivered(forwarder(k)) -1;
                        else
                            state(1) = 2;
                        end
                        t_count(forwarder(k)) = 0;
                    end
                end
                forwarder = find(delivered & (state(1:n) - 1)); % updated list
                % check who can transmit without interference, prefrence
                % will be given to farther away nodes
                if isempty(forwarder) ~= 1
                    tmp = forwarder(length(forwarder));
                    for k = length(forwarder)-1:-1:1
                        if forwarder(k+1)-forwarder(k) > 2*T % interferring range, even in ETX single path, nodes within a larger range should be waiting
                            tmp = [tmp forwarder(k)];
                        end
                    end
                    forwarder = tmp;
                end
                clear tmp
                
                 % generate failed/not state now for all nodes 
                for k = 1:length(state)
                    % time and space
                    if k == 1
                        index_fr_st = state(k) + 2;
                    else
                        index_fr_st = state(k) + 2*state(k-1) - 2;
                        % state(k)- kth node's previous state, 
                        % state(k-1) - k-1th node's current state
                        % 1,1 - index = 1, 1, 2- index - 3 time
                        % 2,1 - index = 2 space, 2,2 - indx = 4
                    end
                    tmp_f(k) = rand <= (1 - fr_st(index_fr_st));
                    state(k) = tmp_f(k) + 1;
                end

                if isempty(forwarder) ~= 1
                    for k = 1:length(forwarder)
                        M = min([forwarder(k)+hop   n+1]);
                        tmp = ((rand <= Prob_tmp(M-forwarder(k))).*tmp_f(M));
                        t_count(forwarder(k)) = t_count(forwarder(k)) + 1;

                        if tmp == 1
                            if forwarder(k) ~= 1 % source always has a packet to send
                                delivered(1,forwarder(k)) = delivered(1,forwarder(k))-1; % packet has been forwarded
                            else
                                state(1) = 2; % Source is alive
                            end
                            t_count(forwarder(k)) = 0; % reset no. of transmissions
                            if M == n + 1
                                nreached = nreached + 1;
                            else
                                delivered(1,M) = delivered(1,M)+1;
                            end
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
        step_time = packet_size/bit_rate + ACK_size/bit_rate;
        capacity_etx_mts(i,height_index,Power_index) = mean(capacity_tmp*packet_size/step_time*1e-6); % in Mbits/s
        capacity_etx_mts_std(i,height_index,Power_index) = std(capacity_tmp*packet_size/step_time*1e-6);
        P_emp_etx_mts(i,height_index,Power_index) = mean(ndelivered/Npackets);%ndelivered/Ntrials/Npackets;
        P_emp_etx_mts_std(i,height_index,Power_index) = std(ndelivered/Npackets);
    end
end