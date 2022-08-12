tol = 0.1;          % Tolerance

E = 200000;         % Young's Modulus
b_s = 500;          % Frame Spacing
sigma_c = 160;      % Characteristic stress level
gamma0 = 1.25;      % Factor of Safety
K1 = 4;             % Plate buckling coefficient
K2 = 0.1;           % Stiffener torsional buckling coefficient
K3 = 0.5;           % Panel buckling coefficient

% Objective (Cost) function that needs to be minimized
% z = f1 * tp + f2 * hs
f1 = 2;
f2 = 1;
f = [f1; f2;];

LB = [5; 20;];      % Lower bounds for design variables
UB = [30; 100;];    % Upper bounds for design variables

x0 = [30; 100;];    % Initial solution
xn = x0;            % Initial solution specification
fval = f'*xn;       % Value of objective function at initial condition

% Iterate until convergence
while true
    
    Dx = 0.05*xn;
    LBn = xn - Dx;
    UBn = xn + Dx;
    
    LBn(LBn < LB) = LB(LBn < LB);
    UBn(UBn > UB) = UB(UBn > UB);
    
    % Linearize constraints    
    A = zeros(3,2);
    b = zeros(3,1);
    
    %----------------------------------------------------------------------
    % 
    % Write the code in the space below to compute the matrix A and b
    %----------------------------------------------------------------------
    t1 = xn(1); h1 = xn(2);   
    A = [-sqrt(K1*E),0; -sqrt(K2*E),sqrt(sigma_c*gamma0); (-K3*E*h1),(-K3*E*t1)];
    b = [sqrt(sigma_c*gamma0)*(-b_s); 0; -(K3*E*h1*t1)-(gamma0*(b_s^2)*sigma_c)]; 
    
    %----------------------------------------------------------------------
    
    % Optimize using linear programming
    [xn1,fvaln,exitflag] = linprog(f,A,b,[],[],LBn,UBn);
    
    if exitflag ~= 1
        error('Optimization failed. Recheck the problem parameters')
    end
    
    xn = xn1;
    
    if abs(fval - fvaln) < tol
        break
    end    
    
    fval = fvaln;    
    
end

fprintf('The optimum solution is tp = %.2f mm and hs = %.2f mm\n',xn(1),xn(2));

seqlinprog_plot
