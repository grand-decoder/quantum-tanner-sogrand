function [l_e, n_guess]=component_decoding(l_in, list_size, sogrand_thr, H_comp, syndrome_comp, IC)
% COMPONENT_DECODING  Performs component decoding using SOGRAND.
% Arguments:
%   l_in           - Input LLRs for the component code (column vector)
%   list_size      - Maximum size of the candidate list for SOGRAND
%   sogrand_thr    - Threshold for the SOGRAND algorithm (between 0 and 1). If the estimated list error probability 
%                    is below this threshold, SOGRAND stops list decoding.
%   H_comp         - Parity-check matrix of the component code
%   syndrome_comp  - Measured syndrome for the component code (column vector) 
%   IC             - Intercept (integer >= 0) for ORBGRAND noise generator. 
%                    If IC = -1, the intercept is automatically determined based on the
%                    current LLRs.
%
% Returns:
%   l_e            - Extrinsic output LLRs after component decoding
%   n_guess        - Number of guesses made by SOGRAND

    % Currently, a fixed high abandonment threshold is used. List decoding is terminated via list_size and
    % sogrand_thr
    Tmax = 2^32;
    [chat_list, s_list, n_guess, curL, pNL, ~] = SOGRAND(l_in, uint8(reshape(H_comp', [], 1)), int32(IC), uint64(list_size), uint64(Tmax), 1-sogrand_thr, false, uint8(syndrome_comp));
    [~, l_e] = list2llr(chat_list, s_list(2, :), curL, pNL, l_in);
end

function [L_APP, L_E] = list2llr(chat, PM, curL, p_notinlist, llr)
% Computes a posterori and extrinsic LLRs from the C SOGRAND output via marginalization.
    chat = chat(:,1:curL);
    PM = PM(:,1:curL);
    assert( size(chat,2) == length(PM), "List size error !")
    p_notinlist = max(p_notinlist, eps);
    
    pp1 = 1./(1+exp(llr));
    pp0 = 1 - pp1;
    
    %
    pp1 = max(pp1, eps); pp1 = min(pp1, 1-eps);
    pp0 = max(pp0, eps); pp0 = min(pp0, 1-eps);
    
    p = exp(-PM);
    
    p1 = sum(chat.*p,2);
    p0 = sum((1-chat).*p,2);
    
    p0 = p0 + p_notinlist * pp0;
    p1 = p1 + p_notinlist * pp1;
    
    L_APP = log(p0) - log(p1);
 
    L_E = L_APP - llr;
end