function [Pf,bet] = safety_index(mu1,sig1,mu2,sig2)

% Write the code to compute the probability of failure assuming there are
% no approximational uncertainties

Pf = integral(@(y) logncdf(y,mu2,sig2).*lognpdf(y,mu1,sig1), 0, Inf);


% Calculate the corresponding safety factor beta for the given probability
% distributions of Q and QL
mean1 = exp(mu1 + ((sig1^2)/2));
mean2 = exp(mu2 + ((sig2^2)/2));

std1 = (exp(sig1^2) - 1)*exp((2*mu1) + (sig1^2));
std2 = (exp(sig2^2) - 1)*exp((2*mu2) + (sig2^2));

bet = (mean2 - mean1)/(sqrt((std1) + (std2)));

end
