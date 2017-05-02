function bits = rlgr(R)
%RLGR computes number of bits to code R with Adaptive Run Length Golomb Rice
%
% R = Nx1 array of integers to be coded
%
% bits = number of bits needed to code R

% Length of input.
N = length(R);

% Constants.
L = 4;
U0 = 3;
D0 = 1;
U1 = 2;
D1 = 1;

% Initialize state.
bits = 0;
k_P = 0;
k_RP = 10*L;

% Preprocess data from signed to unsigned.
U = 2 * R;
neg = (R < 0);
Uneg = -U(neg);
U(neg) = Uneg - 1;

% Process data one sample at a time (time consuming in Matlab).
n = 1;
while n <= N
    
    k = floor(k_P / L);
    k_R = floor(k_RP / L);
    
    u = U(n); % symbol to encode
    
    if k == 0 % no-run mode
        
        % Output GR code for symbol u.
        % bits = bits + gr(u,k_R);
        bits = bits + (floor(u/(2^k_R)) + 1 + k_R);
        
        % Adapt k_R.
        p = floor(u/(2^k_R)); % number of probability halvings
        if p == 0
            k_RP = max(0,k_RP - 2);
        elseif p > 1
            k_RP = k_RP + p + 1;
        end
        
        % Adapt k.
        if u == 0
            k_P = k_P + U0;
        else % u > 0
            k_P = max(0,k_P - D0);
        end
    else % k > 0 % run mode
        
        m = bitshift(1,k); % m = 2^k = expected length of run of zeros
        
        % Parse off up to m symbols,
        % through first non-zero symbol,
        % counting number of zero symbols before it.
        zeroCount = 0;
        while u == 0
            zeroCount = zeroCount + 1;
            if zeroCount >= m || n >= N
                break;
            end
            n = n + 1;
            u = U(n);
        end
        % At this point, either u>0 or (u=0 & (zeroCount>=m | n>=N).
        % That is, either u>0 or (u=0 & zeroCount>=m) or (u=0 & n>=N).
        if zeroCount == m
            % Found a complete run of zeroCount = m zeros.
            % Output a 0.
            bits = bits + 1;
            
            % Adapt k.
            k_P = k_P + U1;
        else % zeroCount < m, and either u>0 or (u=0 and n>=N)
            % Found a partial run of zeroCount < m zeros.
            if u > 0
                % Partial run ended normally with a non-zero symbol u.
                % Output a 1 + length of partial run + GR code for non-zero symbol.
                % bits = bits + 1 + k + gr(u-1,k_R);
                bits = bits + 1 + k + (floor((u-1)/(2^k_R)) + 1 + k_R);
                
                % Adapt k_R.
                p = floor((u-1)/(2^k_R)); % number of probability halvings
                if p == 0
                    k_RP = max(0,k_RP - 2);
                elseif p > 1
                    k_RP = k_RP + p + 1;
                end
                
                % Adapt k.
                k_P = max(0,k_P - D1);
            else % u = 0 and n = N
                % Partial run ended with a zero symbol, at end of sequence.
                % Output a 0.  Leave it to decoder to know the number of symbols needed.
                bits = bits + 1;
            end
        end
    end
    
    n = n + 1;
end

