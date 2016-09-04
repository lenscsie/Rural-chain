% Outage probability calculation using capacity equation
% We are assuming:
% symmetric links - not true always
% random node failures - temporary state if un-responsiveness due to
% protocols like MAC management etc

close all
clear all

% Parameters to set
Distance = 100; % Total range of the wireless chain in km 
%d = 1:1:100; % Distance in km b/w the nodes
% instead of using same d for all heights, fix N_max, such that Prob of
% packet delivery between source and N_max node is >= P_th
N_max = 10;
Prob_th = 0.10; % 1 out of 10 packets is a good threshold for outage, antyhing smaller is too small and could be ignored
T_R = 3; % Number of allowed re-transmissions

Power = [27+12 36+6 27+30 22+60];% point-to-point FCC

% Parameters of path loss model, page 41, Mark
fc = 2.4e3; %MHz
height = [9 15 24 30 45];%[2 9 15 21 24 30 45]; %m
sigma = 8; % Shadowing standard deviation
sigma_awgn_exp = 18; % In fact it is 21 for dBm 1e3 will make it 18
sigma_awgn = 10^(-sigma_awgn_exp)*22e6/2;% (-21) 300K * Boltzmann constant (-11) assumed to yield typical 802.11b ranges %(10^(-170/10)/1000)/2*22e6; % No/2*BW, No = (10^(-x/10)/1000), x is in dBm, BW = 22 MHz
bit_rate = 54e6;
bandwidth = 22e6;
R = bit_rate/bandwidth; %Mbits/s/Hz

packet_size = 1024*8;% in bits
ACK_size = 14*8; % 14 bytes 802.11 standard
Ntrials = 10; % No. of simulation trials
Npackets = 1000; % No. of packets in our simulation trials

% iid node failures
pf = 1e-2; % Probability that a node would fail, same for all nodes, 

% Failures are correlated in Time, Markov model, 
fr = [0.6 1e-2]; % fr(1) - when node was failed in the previous cycle 

% Failures are correlated in Time and Space
fr_st = [0.7 0.5 0.6 1e-2]; % both fail, space fail, time fail, no fail


M_hybrid = 5; % Hybrid chain parameter, positions of redundant peer every mth node
d_peer = 1e-3; % distance between redundant peers in Kms

%n = ceil(Distance./d); % number of links



for Power_index = 1:length(Power)
   for height_index = 1:length(height)
        PtGtGr = Power(Power_index);
        h = height(height_index);

        a =  (1.1*log10(fc) - 0.7)*h - (1.56*log10(fc)-0.8);
        A = 69.55 + 26.16*log10(fc) - 13.82*log10(h) - a;
        B = 44.9 - 6.55*log10(h);
        D = 40.94 + 4.78*(log10(fc))^2 - 18.33*log10(fc);
        alpha = PtGtGr - A + D; % uaing Max EIRP from FCC
        beta = B;

%         for i = 1:max(n) % Very inefficient way :(
%             p_l(i,:) = 0.5*(1 - erf((10*log10((2^R - 1)*sigma_awgn) - alpha + beta.*log10(d*i))/(sqrt(2)*sigma))); % Ps = 1 - P_outage
% 
%         end

        Prob = 0.5*(1 - erf((10*log10((2^R - 1)*sigma_awgn) - alpha + beta.*log10([1:Distance]))/(sqrt(2)*sigma))); % Ps = 1 - P_outage
        flag_mimo = 0;
        
        simple_formula_3_pf_delay_vard % script to calculate delivery probability for simple chain using formula, output P_n

        simple_simulation_opp_effCapacity_vard % output P_emp T_R should be 1.
        
        done = 1
        M_hybrid_old = M_hybrid; % Inorder to run simple chain w/ Markov failures
        M_hybrid = 0;
        hybrid_simulation_opp_effCapacity_MarkovTime_vard 
        
        done = 2
        capacity_m(:,height_index,Power_index) = capacity_hybrid_m(:,height_index,Power_index); % in Mbits/s
        capacity_m_std(:,height_index,Power_index) = capacity_hybrid_m_std(:,height_index,Power_index); % in Mbits/s
        P_emp_m(:,height_index,Power_index) = P_emp_hybrid_m(:,height_index,Power_index);
        P_emp_m_std(:,height_index,Power_index) = P_emp_hybrid_m_std(:,height_index,Power_index);
        
        hybrid_simulation_opp_effCapacity_MarkovTime_Space_vard
        
        done = 3
        capacity_mts(:,height_index,Power_index) = capacity_hybrid_mts(:,height_index,Power_index); % in Mbits/s
        capacity_mts_std(:,height_index,Power_index) = capacity_hybrid_mts_std(:,height_index,Power_index); % in Mbits/s
        P_emp_mts(:,height_index,Power_index) = P_emp_hybrid_mts(:,height_index,Power_index);
        P_emp_mts_std(:,height_index,Power_index) = P_emp_hybrid_mts_std(:,height_index,Power_index);
        M_hybrid = M_hybrid_old;
               
        simple_etx_2_delay_vard % Conventional ETX single path routing

        simple_etx_2_simulation_effCapacity_vard

        %simple_etx_2_simulation_effCapacity_optimistic_vard
        
        simple_etx_2_simulation_effCapacity_MarkovTime_vard
        
        done = 4
        simple_etx_2_simulation_effCapacity_MarkovTime_Space_vard
        
        done = 5
        
        hybrid_formula_3_pf_delay_vard

        hybrid_simulation_opp_effCapacity_vard
        
        done = 6

        %hybrid_formula_3_pf_delay_MarkovTime_3
        hybrid_simulation_opp_effCapacity_MarkovTime_vard 
        
        done = 7
        hybrid_simulation_opp_effCapacity_MarkovTime_Space_vard
        
        done = 8
        filename = strcat('results_g_M',num2str(M_hybrid),'_TR',num2str(T_R),'_',num2str(sigma_awgn_exp),'_h',num2str(h),'_P',num2str(PtGtGr));
        save(filename);
        done = 10
    end
end
% P_n_hybrid_m0(:,1) = P_n_hybrid(:,height_index,Power_index);
% P_emp_hybrid_m0(:,1) = P_emp_hybrid(:,height_index,Power_index);
% P_emp_hybrid_m0_std(:,1) = P_emp_hybrid_std(:,height_index,Power_index);
% 
% save results_g_M0_TR3_h9_P45 P_n_hybrid_m0 P_emp_hybrid_m0 P_emp_hybrid_m0_std
% figure
% title('Reliability')
% plot([1:N_max],P_n,[1:N_max],P_emp,[1:N_max],P_n_etx,[1:N_max],P_emp_etx)
% hold on
% plot(P_emp_etx_m)
% plot(P_emp_etx_mts)
% plot(P_emp_m)
% plot(P_emp_mts)
% plot(P_emp_etx_op)
% plot(P_n_hybrid)
% plot(P_emp_hybrid_m)
% plot(P_emp_hybrid_mts)
% 
% % xlabel('d (km)')
% ylabel('P_n')
% % legend('Opportunistic routing (OR)','OR Simulation','ETX', 'ETX Simulation')
% 
% figure
% title('Effective data rate')
% plot(d(1:length(capacity)),capacity,d(1:length(capacity_etx)),capacity_etx,d(1:length(capacity_etx_op)),capacity_etx_op)
% xlabel('d (km)')
% ylabel('Data rate (Mbits/s)')
% legend('Opportunistic routing (OR)','ETX','ETX optimistic' )
% 
% figure
% title('Average end-to-end delay')
% plot(d(1:length(avg_delay)),avg_delay,d(1:length(avg_delay_etx)),avg_delay_etx)
% xlabel('d (km)')
% ylabel('Average end-to-end delay (sec)')
% legend('Opportunistic routing (OR)','ETX')

% figure
% title('Reliability')
% plot(d(1:length(P_n)),P_n,d(1:length(P_emp)),P_emp,d(1:length(P_n_hybrid)),P_n_hybrid,d(1:length(P_emp_hybrid)),P_emp_hybrid)
% %plot(d(1:length(P_n)),P_n,d(1:length(P_n_hybrid)),P_n_hybrid)
% plot(d(1:length(P_n_hybrid)),P_n_hybrid,d(1:length(P_emp_hybrid)),P_emp_hybrid)
% xlabel('d (km)')
% ylabel('P_n')
% legend('simple theoretical','simple simulation',strcat('Hybrid m = ',num2str(M_hybrid),' theoretical'),strcat('Hybrid m = ',num2str(M_hybrid),' simulation'))
% %legend(strcat('Hybrid m = ',num2str(M_hybrid),' theoretical'),strcat('Hybrid m = ',num2str(M_hybrid),' simulation'))
% 
% 
% figure
% title('Effective data rate')
% plot(d(1:length(capacity)),capacity,d(1:length(capacity_hybrid)),capacity_hybrid)
% xlabel('d (km)')
% ylabel('Data rate (Mbits/s)')
% legend('Simple',strcat('Hybrid m = ',num2str(M_hybrid)))
% 
% figure
% title('Average end-to-end delay')
% plot(d(1:length(avg_delay)),avg_delay,d(1:length(avg_delay_hybrid)),avg_delay_hybrid)
% xlabel('d (km)')
% ylabel('Average end-to-end delay (sec)')
% legend('Simple',strcat('Hybrid m = ',num2str(M_hybrid)))
