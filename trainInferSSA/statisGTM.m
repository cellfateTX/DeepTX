function statis = statisGTM(param,k)
%% This code computes some statistics of the GTM.
% Inputs:
%    param: A structure contains model parameters 
%    param.kon and param.ron: OFF dwell time distribution f_off(t) = ron^(kon) * t^(kon-1) * e^(-ron * t) / gamma(kon)
%    param.koff and param.roff: ON dwell time distribution f_on(t) = roff^(koff) * t^(koff-1) * e^(-roff * t) / gamma(koff)
%    param.mu: Transcriptional rate
%    param.delta: Degradation rate.
%
% Output: statis: [mean,noise,sk,kt,bc,ff].

% parameters setting
kon = param.kon;
ron = param.ron;
koff = param.koff;
roff = param.roff;
mu = param.mu;
delta = param.delta;

Laplace_s = 0:k;
Lfon = (ron./(Laplace_s + ron)).^kon; % Laplace fon(x)
Lfoff = (roff./(Laplace_s + roff)).^koff;% Laplace foff(x)
mean_tau_off = kon/ron;
mean_tau_on = koff/roff;
c = 1/(mean_tau_off + mean_tau_on);
LFon = [mean_tau_off (1-Lfon(2:end))./Laplace_s(2:end)]; % Laplace Fon(x)
LFoff = [mean_tau_on (1-Lfoff(2:end))./Laplace_s(2:end)]; % Laplace Foff(x)

Ck = zeros(1,k);
Ck(1) = (1 - Lfoff(2))/(mean_tau_off + mean_tau_on);
for iter = 2:k
    i = 1:iter;
    Ck_coef = [(mu.^i(1:end-1)) .* Lfon(iter-i(1:end-1)+1) .* Ck(iter-i(1:end-1)) ./ ((1 - Lfoff(iter-i(1:end-1)+1) .* Lfon(iter-i(1:end-1)+1)) .* gamma(i(1:end-1)+1)), c*mu^(iter-1)/gamma(iter+1)];
    Ck_sum = zeros(1,iter);
    for iter_i = 1:iter
        for iter_j = 0:iter_i
            Ck_sum(iter_i) = Ck_sum(iter_i) + matnchoosek(iter_i,iter_j) * (-1)^(iter_i-iter_j) * Lfoff(iter-iter_j+1);
        end
    end
    Ck(iter) = sum(Ck_coef .* Ck_sum);
end

bk = zeros(1,k);
for iter = 1:k
    i = 1:iter;
    bk_coef = [c*mu^(iter)/gamma(iter+1), (mu / iter) * (mu.^(iter-i(1:end-1))) .* Lfon(i(1:end-1)+1) .* Ck(i(1:end-1)) ./ (gamma(iter-i(1:end-1)) .* (1 - Lfon(i(1:end-1)+1) .* Lfoff(i(1:end-1)+1)))];
    bk_sum = zeros(1,iter);
    for iter_i = 0:(iter - 1)
        for iter_j = 0:(iter - iter_i -1)
            bk_sum(iter_i+1) = bk_sum(iter_i+1) + matnchoosek(iter-iter_i-1,iter_j) * (-1)^(iter-iter_i-iter_j-1) * LFoff(iter-iter_j);
        end
    end
    bk(iter) = sum(bk_coef .* bk_sum);
end

b = bk;
m = b(1);
v = 2*b(2)+b(1)-b(1)^2;
cv2 = v/(m^2);
fano = v/m;
sk = (6*b(3)+6*b(2)+b(1)-3*b(1)*(2*b(2)+b(1))+2*b(1)^3)/(2*b(2)+b(1)-b(1)^2)^1.5 + 1;
kt = (24*b(4)+36*b(3)+14*b(2)+b(1)-4*b(1).*(6*b(3)+6*b(2)+b(1))+...
    6*b(1)^2*(2*b(2)+b(1))-3*b(1)^4)/(2*b(2)+b(1)-b(1)^2)^2;
bc = ((sk - 1)^2 + 1)/kt;
statis = [m,cv2,fano,sk,kt,bc];
end

function c = matnchoosek(n,k) 
        g = gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1);
        c = floor(exp(g)+.5);
    end
