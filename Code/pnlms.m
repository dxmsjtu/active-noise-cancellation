function [y] = pnlms(echo_file, far_file, out_file, taps, myu, dump_flag)
% [1]Proportionate Normalized Least-Mean-Squares Adaptation in Echo Cancelers 
    [echo, sr_echo] = audioread(echo_file);    [far, sr_far] = audioread(far_file);
    if(sr_echo~=sr_far)         error('echo signal sampling rate must be equal to far signal');    end
    % para
    min_len = min(length(echo),length(far));    out_len = min_len-taps;
    out = zeros(out_len, 1);    w = zeros(taps, 1);    sigma = 0.01;    rou = 0.01;
    for i=1:out_len
        x = fliplr(far(i:i+taps-1));
        d = echo(i);
        out(i) = d - sum(x.*w);
        e = out(i);        
        l_inf = max(abs(w));                     % (9) of [1] 
        l_prime_inf = max(sigma, l_inf);  % (10) of [1] 
        g = max(rou*l_prime_inf, abs(w));% (11) of [1] 
        g_norm = sum(g)/taps + eps;      % (12) of [1] 
        norm = x' * x + eps;                    % (13) of [1]  
        w = w + g ./ g_norm .* myu .* e .* x ./ norm;% (14) of [1] 
    end
    if(dump_flag)
        audiowrite(out_file, out, sr_echo); % figure; plot(far); figure; plot(out,'r');
    end
    y = erle(out, echo(1:out_len));
end
