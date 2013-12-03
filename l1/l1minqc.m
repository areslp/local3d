function xp = l1minqc(K, f, y, epsilon, x0, lytol, mu, newtonmaxiter, maxiter, itertol) 
% xp = l1minqc(K, f, y, epsilon, x0, lytol, mu, newtonmaxiter, maxiter, itertol) 
%
% Solve quadratically constrained l1 minimization:
% min ||K * x - f||_1   s.t.  ||y - x||_2 <= \epsilon
%
% Reformulate as the second-order cone program
% min_{x,u}  sum(u)   s.t.    Kx - u - f<= 0,
%                            -Kx - u + f<= 0,
%      1/2(||y - x||^2 - \epsilon^2) <= 0
% and use a log barrier algorithm.
%  
%
% x0 - Nx1 vector, initial point.
%
% K - Matrix K which is m-by-n. m >= n. Three options:
%     dense matrix - either \ or blendenpik is used to solve linear
%                    equations (for search direction).
%     sparse matrix - use LSQR to solve linear equations using a 
%                     preconditionning.
%     function handle - represents Kx and K'x. 
%                       K * x = K(x, 'notransp')
%                       K' * x = K(x, 'transp')
%                       Use unpreconditioned LSQR. 
%
%
% f - mx1 right hand side.
%
% y - nx1 base solution.
%
% epsilon - scalar - max distance from base solution.
%
% lytol - The log barrier algorithm terminates when the duality gap <= lytol.
%         Also, the number of log barrier iterations is completely
%         determined by lytol.
%         Default = 1e-3.
%
% mu - Factor by which to increase the barrier constant at each iteration.
%      Default = 10.
%
% newtonmaxiter - Max iterations for Newton method.
%                 Default = 20.
%
% maxiter - Maximum number of iterations for LSQR.
%           Ignored if K is dense.
%           Default = 100.
%
% itertol - Tolerance for LSQR for sparse K.
%           Ignored if K is dense.
%           Default = 1e-8.
%
% Algorithm is described in Section 5 of:
% "L1-sparse reconstruction of sharp point set surfaces"
% http://www.cs.tau.ac.il/~haima/l1sparse-tog-final.pdf
%
% Written by: Haim Avron, based on l1-Magic
% Email: haim.avron@gmail.com
% Created: April 2009
%

%% Complete parameters
if (nargin < 6)
    lytol = 1e-3; 
end
if (nargin < 7)
    mu = 10; 
end
if (nargin < 8) 
    newtonmaxiter = 20; 
end
if (nargin < 9) 
    maxiter = 100; 
end
if (nargin < 10)
    itertol = 1e-8; 
end

tstart = cputime;

funchand = isa(K,'function_handle');

newtontol = lytol;
N = length(x0);

x = x0;
if (funchand)
    d = abs(K(x0, 'notransp') - f);
else
    d = abs(K * x0 - f);
end
u = (0.95)*(d) + (0.10)*max(d);

disp(sprintf('Original l1 norm = %.3f, original functional = %.3f', sum(d), sum(u)));

% choose initial value of tau so that the duality gap after the first
% step will be about the origial norm
tau = (2 * length(d) + 1) / sum(d);

lyiter = ceil((log(2 * length(d) + 1)-log(lytol)-log(tau))/log(mu));
disp(sprintf('Number of log barrier iterations = %d\n', lyiter));

totaliter = 0;

for ii = 1:lyiter

  [xp, up, ntiter] = l1minqc_newton(x, u, K, f, y, epsilon, tau, newtontol, newtonmaxiter, itertol, maxiter);
  totaliter = totaliter + ntiter;
  
  if (funchand)
      res = abs(K(xp, 'notransp') - f);
  else
      res = abs(K * xp - f);
  end
  disp(sprintf('\nLog barrier iter = %d, l1 = %.3f, functional = %8.3f, tau = %8.3e, total newton iter = %d\n',...
       ii, sum(res), sum(up), tau, totaliter));
  
  x = xp;
  u = up;
 
  tau = mu*tau;
end

solvetime = cputime - tstart;
disp(sprintf('\n\nTotal solution time (CPU time) is %.2e seconds.',  solvetime));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Newton method
function [xp, up, niter] = l1minqc_newton(x0, u0, K, f, y, epsilon, tau, newtontol, newtonmaxiter, itertol, maxiter)

%% Determine external libraries
use_blendenpik = ~isempty(which('blendenpik'));
use_sptrisolve = ~isempty(which('sptrisolve'));

if (use_blendenpik)
    blendenpik_params.tol = itertol;
end

%% Determine matrix type
if isa(K, 'function_handle');
    Ktype = 2;
    Kapply = @(x) K(x, 'notransp');
    Ktapply = @(x) K(x, 'transp');
else
    Kt = K';
    Kapply = @(x) K * x;
    Ktapply = @(x) Kt * x;
    if (issparse(K))
        Ktype = 1;
    else
        Ktype = 0;
    end
end

% line search parameters
alpha = 0.01;
beta = 0.5;

N = length(x0);
M = length(f);

% initial point
x = x0;
u = u0;
r = x - y;
fu1 = Kapply(x) - u - f;
fu2 = -Kapply(x) - u + f;
% if (Ktype == 2)
%     fu1 = K(x, 'notransp') - u - f;
%     fu2 = -K(x, 'notransp') - u + f;
% else
%     fu1 = K * x - u - f;
%     fu2 = - K * x - u + f;
% end
fe = 1/2*(r'*r - epsilon^2);
fc = sum(u) - (1/tau)*(sum(log(-fu1)) + sum(log(-fu2)) + log(-fe));


perm = [];

niter = 0;
done = 0;
while (~done)

    ntgz = 1./fu1 - 1./fu2;
    ntgu = -tau - 1./fu1 - 1./fu2;

    sig11 = 1./fu1.^2 + 1./fu2.^2;
    sig12 = -1./fu1.^2 + 1./fu2.^2;
    %% The following are the 'true' formulas, but to avoid overflow other
    %% formulas are used.
    %sigx = sig11 - sig12.^2./sig11;
    %ssigx = sqrt(sigx);
    ssigx = 2 ./ (abs(fu1) .* abs(fu2) .* sqrt(sig11));
    SSigx = spdiags(ssigx, 0, M, M);
% 
%     if (Ktype == 2)
%         gradf = -(1/tau)*[K(ntgz, 'transp') + 1/fe*r; ntgu];
%     else
%         gradf = -(1/tau)*[Kt * ntgz + 1/fe*r; ntgu];
%     end
    gradf = -(1/tau)*[Ktapply(ntgz) + 1/fe*r; ntgu];
    w = [(ntgz - sig12./sig11.*ntgu) ./ ssigx;
        zeros(N, 1);
        1];
    
    % Now we need to solve min(norm(Ktilde*dx - w, 2)) 
    tstart = cputime;
    switch(Ktype) 
        case 0
            Ktilde = [SSigx * K; (1/sqrt(-fe)) * eye(N, N); (1/fe) * r'];
            if (max(ssigx) / (-fe) < 1e3)
                if (use_blendenpik)
                    dx = blendenpik(Ktilde, w, blendenpik_params);
                    solvemethod = 'Blendenpik';
                else
                    dx = Ktilde \ w;
                    solvemethod = 'Direct (MATLAB backslash)';
                end
            else
                % In this case Ktilde should be well-conditioned to begin
                % with. We can use an iterative method without worry.
                solvemethod = 'Iterative';
                [dx, flag] = lsqr(Ktilde, w, itertol, maxiter);
                if (flag ~= 0)
                    disp('WARNING: LSQR failed to converege!!!');
                end
            end
            
        case 1
                
            Ktilde1 = [SSigx * K; (1/sqrt(-fe)) * speye(N, N)];
            Ktilde = [Ktilde1; (1/fe) * r'];

            if (max(ssigx) / (-fe) < 1e3)
                solvemethod = 'Semi-direct Cholesky';
                if (~isempty(perm))
                    Ktilde1P = Ktilde1(:, perm);
                    clear L;
                    [L, p] = chol2(Ktilde1P' * Ktilde1P);
                else
                    [L, p, perm] = chol2(Ktilde1' * Ktilde1);
                    if (p > 0)
                        solvemethod = 'Semi-direct QR';
                        opts.Q = 'discard';
                        opts.econ = N;
                        clear L;
                        [Qdummy, L, P] = spqr(Ktilde1, opts);
                        [perm, dummy] = find(P);
                        p = 0;
                    end
                end

                if (p > 0)
                    solvemethod = 'Semi-direct QR';
                    clear L;
                    L = spqr(Ktilde1P, 0);
                end
                % Using sptrisolve explicitly to solve the traingular systems
                % is much faster! 
                if (use_sptrisolve)
                    precond_fn = @(b, transp) sptrisolve(L, b, 'upper', transp); 
                else
                    precond_fn = @(b, transp) precond_by_L(L, b, transp);
                end
                [dx1, flag] = lsqr(Ktilde(:, perm), w, itertol, maxiter, precond_fn);
                if (flag ~= 0)
                    disp('WARNING: LSQR failed to converege!!!');
                end
                dx(perm, 1) = dx1;
            else
                % In this case Ktilde should be well-conditioned to begin
                % with. We can use an iterative method without worry.
                solvemethod = 'Iterative';
                [dx, flag] = lsqr(Ktilde, w, itertol, maxiter);
                if (flag ~= 0)
                    disp('WARNING: LSQR failed to converege!!!');
                end
            end

        case 2
            solvemethod = 'Iterative';
            Ktilde = @(w, transp) large_scale_Ktilde(K, N, fe, r, SSigx, w, transp);
            [dx, flag, relres, iter] = lsqr(Ktilde, w, itertol, maxiter);
            if (flag ~= 0)
                disp('WARNING: LSQR failed to converege!!!');
            end
    end
    solvetime = cputime - tstart;
    disp(sprintf('\t\tEquation solve time (CPU time): %.2e\tMethod: %s', solvetime, solvemethod));

    Kdx = Kapply(dx);
    du = (ntgu - sig12 .* Kdx) ./ sig11;

    % minimum step size that stays in the interior
    s = 1;
    xp = x + s*dx;  up = u + s*du;  rp = r + s * dx;
    coneiter = 0;
    while ( (max(abs(Kapply(xp) - f)-up) > 0) || (rp'*rp > epsilon^2) )
        s = beta*s;
        xp = x + s*dx;  up = u + s*du;  rp = r + s*dx;
        coneiter = coneiter + 1;
        if (coneiter > 32)
            disp('Stuck on cone iterations, returning previous iterate.');
            xp = x;  up = u;
            return
        end
    end

    % backtracking line search
    %% TODO: func handle
    fu1p = Kapply(xp) - f - up;  fu2p = f - Kapply(xp) - up;  
    fep = 1/2*(rp'*rp - epsilon^2);
    fcp = sum(up) - (1/tau)*(sum(log(-fu1p)) + sum(log(-fu2p)) + log(-fep));
    fclin = fc + alpha*s*(gradf'*[dx; du]);
    backiter = 0;
    while (fcp > fclin)
        s = beta*s;
        xp = x + s*dx;  up = u + s*du;  rp = r + s*dx;
        fu1p = Kapply(xp) - f - up;  fu2p = f - Kapply(xp) - up;  fep = 1/2*(rp'*rp - epsilon^2);
        fcp = sum(up) - (1/tau)*(sum(log(-fu1p)) + sum(log(-fu2p)) + log(-fep));
        fclin = fc + alpha*s*(gradf'*[dx; du]);
        backiter = backiter + 1;
        if (backiter > 32)
            disp('Stuck on backtracking line search, returning previous iterate.');
            xp = x;  up = u;
            return
        end
    end

    % set up for next iteration
    x = xp; u = up;  r = rp;
    fu1 = fu1p;  fu2 = fu2p;  fe = fep;  fc = fcp;

    lambda2 = -(gradf'*[dx; du]);
    stepsize = s*norm([dx; du]);
    niter = niter + 1;
    done = (lambda2/2 < newtontol) | (niter >= newtonmaxiter);

    disp(sprintf('Newton iter = %d, Functional = %.2e, Newton decrement = %.2e, Stepsize = %.2e, Cone iterations = %d, Backtrack iterations = %d', ...
        niter, fc, lambda2/2, stepsize, coneiter, backiter));

end
    
%%
function x = precond_by_L(L, b, transp)

if (strcmp(transp, 'transp'))
    x = L' \ b;
else
    x = L \ b;
end

%%
function x = large_scale_Ktilde(K, N, fe, r, SSigx, w, transp)

if (strcmp(transp, 'notransp'))
    x = [SSigx * K(w, 'notransp');
         (1/sqrt(-fe)) * w;
         (1/fe) * r' * w];
else
    M = size(SSigx, 1);
    x = K(SSigx * w(1:M), 'transp') + (1/sqrt(-fe)) * w(M+1:M+N) + (1/fe) * w(end) * r;
end


            
                   
