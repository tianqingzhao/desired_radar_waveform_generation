function s = signal_bpsk_ori(M)

    s = zeros(M, 2^M);
    for i = 1:M
        flag_num = 1;
        flag = 1;
        for j = 1:2^M
           s(i, j) = flag;
           if (flag_num < 2^(M-i))
               flag_num = flag_num + 1;
           else
               flag_num = 1;
               flag = -flag;
           end
        end
    end

end