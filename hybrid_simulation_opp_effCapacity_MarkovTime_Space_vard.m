% Calculates empirical delivery probability for hybrid chain
% Effective capacity in saturation condition
% Not truly accurate: we assume, pkts received at the redundant peer is
% received by the primary node - June 28th, 2012

% fr_st = [0.7 0.5 0.6 1e-2]; % both fail, space fail, time fail, no fail

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
         % For redundant peers at mth position
        Prob_redundant = 0.5*(1 - erf((10*log10((2^R - 1)*sigma_awgn) - alpha + beta.*log10(((d.*[1:T]).^2 + d_peer^2).^(0.5)))/(sqrt(2)*sigma))); % Ps = 1 - P_outage
        Prob_peer = 0.5*(1 - erf((10*log10((2^R - 1)*sigma_awgn) - alpha + beta.*log10(d_peer))/(sqrt(2)*sigma)));

        P_r = Prob_redundant([1:T]); % For redundant peers

        capacity_tmp = zeros(1,Ntrials); % reached packets/time taken
        for j = 1:Ntrials
            delivered = zeros(1,n); % status of nodes having packet or not, 0 if forwarded
            delivered_r = zeros(1,n); % status of redundant peers according to the location
            state = 2*ones(1,n+1); % node is failed or not 1 - failed, 2- not failed
            state_r = 2*ones(1,n);
            delivered(1,1) = 1;
            
            t_count = zeros(1,n); % No. of trans from each node
            t_count_r = zeros(1,n);
            nreached = 0; %no. of packets reaching the destination
            ndropped = 0; %no. of packets dropped
            steps = 0;
            txcount = 1; % Used only for trouble shooting
            
            while(nreached+ndropped < Npackets)

                % forwarder list for this round
                forwarder = find(delivered | delivered_r);
                % re-transmission consideration
                for k = 1:length(forwarder)
                    if (delivered(forwarder(k)) >= 1 && t_count(forwarder(k)) >= T_R) ||...
                            (delivered_r(forwarder(k)) >= 1 && t_count_r(forwarder(k)) >= T_R) % drop the packet
                        ndropped = ndropped + 1;

                        if delivered(1,forwarder(k)) >= 1
                            if forwarder(k) ~= 1
                                delivered(1,forwarder(k)) = delivered(1,forwarder(k))-1;
                            else
                                txcount = txcount+1;
                                state(1) = 2;
                            end
                            t_count(forwarder(k)) = 0;
                        else
                            delivered_r(1,forwarder(k)) = delivered_r(1,forwarder(k)) - 1;
                            t_count_r(forwarder(k)) = 0;
                        end


                    end
                end
                forwarder = find((delivered & (state(1:n) - 1)) | (delivered_r & (state_r - 1))); % updated list
                % check who can transmit without interference, preferrence
                % will be given to farther away nodes
                
                if isempty(forwarder) ~= 1
                    tmp = forwarder(length(forwarder));
                    for k = length(forwarder)-1:-1:1
                        if tmp(length(tmp))-forwarder(k) > 2*T % interferring range
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
                for k = 2:length(state_r)
                    % time and space
                    
                    index_fr_st_r = state_r(k) + 2*state(k-1) - 2;
                    tmp_f_r(k) = rand <= (1 - fr_st(index_fr_st_r));
                    state_r(k) = tmp_f_r(k) + 1;
                end
                

                if isempty(forwarder) ~= 1
                    for k = 1:length(forwarder)
                        tmp_peer = 0;
                        M = min([T   n+1-forwarder(k)]);

                        M2 = min([T   n-forwarder(k)]); % Destination does not have a redundant peer

                        % Determining how many redundant peer we have from
                        % forwarder(k) to fowarder(k) + M
                        index = find(rem(forwarder(k)-1+[1:M2],M_hybrid)== 0);
                        
                        if delivered(forwarder(k)) >= 1 % Primary chain will transmit

                            tmp_r = zeros(1,M2);
                            if isempty(index) ~= 1
                                tmp_r(index) = ((rand(1,length(index)) <= P_r(index)).*tmp_f_r(forwarder(k)+index));

                            end
                            tmp = ((rand(1,M) <= P(1:M)).*tmp_f(forwarder(k) + [1:M])) ;
                        else % Redundant node will transmit
                            tmp_r = zeros(1,M2);
                            if isempty(index) ~= 1
                                tmp_r(index) = ((rand(1,length(index)) <= P(index)).*tmp_f_r(forwarder(k) + index));
                            end
                            tmp = ((rand(1,M) <= P_r(1:M)).*tmp_f(forwarder(k) + [1:M])) ;
                            tmp_peer = ((rand <= Prob_peer).*tmp_f(forwarder(k))) ;

                        end

                        if delivered(forwarder(k)) >= 1
                            t_count(forwarder(k)) = t_count(forwarder(k)) + 1;
                        else
                            t_count_r(forwarder(k)) = t_count_r(forwarder(k)) + 1;
                        end


                        I_del = find(tmp);
                        I_del_r = find(tmp_r);
                        if isempty(I_del) ~= 1 || isempty(I_del_r) ~= 1 || tmp_peer

                            if delivered(forwarder(k)) >= 1
                                if forwarder(k) ~= 1 % source always has a packet to send
                                    delivered(1,forwarder(k)) = delivered(1,forwarder(k))-1; % packet has been forwarded
                                else
                                    txcount = txcount+1;
                                    state(1) = 2; % Source is alive
                                end
                                t_count(forwarder(k)) = 0; % reset no. of transmissions
                            else
                                delivered_r(1,forwarder(k)) = delivered_r(1,forwarder(k))-1; % packet has been forwarded
                                t_count_r(forwarder(k)) = 0; % reset no. of transmissions
                            end

                            if isempty(I_del) ~= 1
                                if forwarder(k) + I_del(length(I_del)) == n + 1
                                    nreached = nreached + 1;
                                else
                                    if isempty(I_del_r) ~= 1
                                        if I_del(length(I_del)) >= I_del_r(length(I_del_r))
                                            delivered(1,forwarder(k) + I_del(length(I_del))) = delivered(1,forwarder(k) + I_del(length(I_del)))+1;
                                        else
                                            delivered_r(1,forwarder(k) + I_del_r(length(I_del_r))) = delivered_r(1,forwarder(k) + I_del_r(length(I_del_r)))+1;
                                        end
                                    else
                                        delivered(1,forwarder(k) + I_del(length(I_del))) = delivered(1,forwarder(k) + I_del(length(I_del)))+1;
                                    end
                                end
                            elseif isempty(I_del_r) ~= 1
                                delivered_r(1,forwarder(k) + I_del_r(length(I_del_r))) = delivered_r(1,forwarder(k) + I_del_r(length(I_del_r)))+1;
                            elseif tmp_peer == 1
                                delivered(1,forwarder(k)) = delivered(1,forwarder(k))+1;
                            end
                        end
                    end
                end
                if nreached >= 1 % just to cut of the time taken for first packet.
                    steps = steps+1;
                end

            end
            %txcount
            %sum(delivered) + sum(delivered_r) + nreached + ndropped
            ndelivered(j) = nreached;%ndelivered = ndelivered + nreached;
            if nreached == 0
                capacity_tmp(j) = 0;
            else
                capacity_tmp(j) = nreached/steps;
            end
        end

        % To calculate capacity in terms of bit/s
        % time taken in each step
        if M_hybrid ~= 0
            step_time = packet_size/bit_rate + (T + ceil(T/M_hybrid) )*ACK_size/bit_rate;
        else
            step_time = packet_size/bit_rate + (T )*ACK_size/bit_rate;
        end
        capacity_hybrid_mts(i,height_index,Power_index) = mean(capacity_tmp*packet_size/step_time*1e-6); % in Mbits/s
        capacity_hybrid_mts_std(i,height_index,Power_index) = std(capacity_tmp*packet_size/step_time*1e-6);
        P_emp_hybrid_mts(i,height_index,Power_index) = mean(ndelivered/Npackets);%ndelivered/Ntrials/Npackets;
        P_emp_hybrid_mts_std(i,height_index,Power_index) = std(ndelivered/Npackets);
    end
end

