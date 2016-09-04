% Traditional ETX based routing
% We are not calculating ETX, rather using the probabilities of links,
% it is logically the same

% Prob = 0.5*(1 - erf((10*log10((2^R - 1)*sigma_awgn) - alpha + beta.*log10([1:Distance]))/(sqrt(2)*sigma))); % Ps = 1 - P_outage
I_p = find(Prob >= Prob_th);


for i = N_max:-1:1% 1:length(d);
    d = I_p(length(I_p))/i; % Internodal distance
    n = ceil(Distance/d); % No of nodes
    if flag_mimo == 1
        Prob_tmp = 0.5*(1 - erf((10*log10((2^(R/Nmin) - 1)*sigma_awgn*Nmin/Nrecv) - alpha + beta.*log10(d*[1:n]))/(sqrt(2)*sigma)));%Prob(I); % Probability of reception at nodes +1,..., +n
    else
        Prob_tmp = 0.5*(1 - erf((10*log10((2^R - 1)*sigma_awgn) - alpha + beta.*log10(d*[1:n]))/(sqrt(2)*sigma)));%Prob(I); % Probability of reception at nodes +1,..., +n
    end
   
    Prob_tmp = Prob_tmp*(1-pf); % Random node failures

    % Re-transmissions
    Prob_link = Prob_tmp;
    for j = 1:T_R-1
        Prob_link = Prob_link + Prob_tmp.*(1-Prob_tmp).^j;
    end
    I = find(Prob_link > 1e-3);
    if isempty(I) ~= 1
        T = I(length(I));
        Prob_link = Prob_link(I);

        hop1 = [1 1];
        preference = [1];
        for j = 2:T
            if Prob_link(j) > prod(Prob_link(hop1))
                hop1 = j;
                preference = [j preference];
            end
            if j < T
                hop1 = [hop1 1];
            end
            % No need to fix the tail here, as we are only interested in the
            % max. hop with highest probability i.e. max(hop1);

        end

        hop = max(hop1);
        K = floor(n/hop);
        tail = rem(n,hop);
        step_time = packet_size/bit_rate + ACK_size/bit_rate;
        step = 0;
        if tail > 0
            P_n_etx(i,height_index,Power_index) = Prob_link(hop)^K*Prob_link(tail);
            step2 = 0;
            for j = 0:T_R-1
                step = step + Prob_tmp(hop)*(1-Prob_tmp(hop)).^j*(j+1)*step_time;
                step2 = step2 + Prob_tmp(tail)*(1-Prob_tmp(tail)).^j*(j+1)*step_time;
            end
            avg_delay_etx(i,height_index,Power_index) = (K)*step+step2;%_time;
        else
            P_n_etx(i,height_index,Power_index) = Prob_link(hop)^K;
            for j = 0:T_R-1
                step = step + Prob_tmp(hop).*(1-Prob_tmp(hop)).^j*(j+1)*step_time;
            end
            avg_delay_etx(i,height_index,Power_index) = (K+1)*step;%_time;
        end
        % tail < hop, if hop gives the best probability, for tail units, tail will give the best probability
        % hop = 6, tail = 4, if Prob_link(6) is the best option, Prob_link(4)
        % will be the best option for tail. Prob_link is a monotonically
        % decreasing function of distance.
        
        n_etx(i) = K+1;
    end
end

