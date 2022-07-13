% Set plotting parameters
set(0,'DefaultAxesFontSize',14)
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesFontName','Times New Roman')

% Parameters describing the PDF of load effect Q
% Note that these are not the mean and standard deviation of load effect Q
mu1 = 0; sig1 = 0.4;

% Parameters describing the PDF of the limit value of the load effect QL
% Note that these are not the mean and standard deviation of load effect QL
mu2 = 1; sig2 = 0.15;

% Probability of Exceedance
p = 0.001:0.001:0.3;

% Array initialization
Qc = zeros(size(p));
QLc = zeros(size(p));

for i = 1:numel(p)
    % Characteristic value of Q corresponding to a probability of
    % nonexceedance of (1 - p(i))
    Qc(i) = logninv(1-p(i),mu1,sig1);
    
    % Characteristic value of QL corresponding to a probability of
    % exceedance of p(i)
    QLc(i) = logninv(p(i),mu2,sig2);
end

plot(p,Qc,p,QLc)
box off
title('Characteristic values')
xlabel('Probability of exceedance $p_e$','interpreter','latex')
ylabel('$Q_c$ and $Q_{L,c}$','interpreter','latex')
legend('$Q_c$','$Q_{L,c}$','interpreter','latex')
legend boxoff

% Calculate the probability of exceedance corresponding to the intersection
% of the two curves above. Make sure to use the FZERO function.

intpfun = @(x) logninv(1-x,mu1,sig1)-logninv(x,mu2,sig2); 
p_int = fzero(intpfun, [0.01,0.5]);

fprintf('Probability of exceedance corresponding to the intersection: %.6f\n',p_int);

% Find the minimum probability of exceedance p such that the system is able
% to accommodate partial safety factors provided below. Make sure to use
% the FZERO function
gammaS1 = 1.1;  % Safety load factor
gammaS2 = 1.2;  % Servicability load factor
gammaQ  = 1.1;  % Approximational load factor
gammaQL = 1.1;  % Approximational usage factor

fun = @(x)(gammaS1*gammaS2*gammaQ*logninv(1-x,mu1,sig1)-logninv(x,mu2,sig2)/gammaQL);
p_psf = fzero(fun,[0.15,0.2]);


fprintf('Probability of exceedance that can accommodate the specified partial safety factors: %.6f\n',p_psf);

% Calculate the area of the overlap of the load and limit state PDFs and 
% compare it with the failure probability in the previous problem. Are the 
% two same?

pQ = @(x) lognpdf(x,mu1,sig1);                          % PDF of load Q
pL = @(x) lognpdf(x,mu2,sig2);                          % PDF of load limit value QL
pT = @(x) lognpdf(x,mu2,sig2)-lognpdf(x,mu1,sig1);
Qstar = fzero(pT,[1,5]);                                                % Intersection point of both PDF curves (Use FZERO function to find this intersection point
                                                        % of PDF curves)
AR = logncdf(Qstar,mu2,sig2) + 1-logncdf(Qstar,mu1,sig1);                                                   % Area of overlap

[Pf,bet] = safety_index(mu1,sig1,mu2,sig2);             % Calculate failure probability using function developed in previous problem

fprintf('Area of overlap: %.6f\n',AR);
fprintf('Probability of Failure: %.6f\n',Pf);
