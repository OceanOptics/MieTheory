% FASTMIE Scattering and absorption of light by homogenous spheres
% [S1,S2] = FASTMIE(X,NREL,ANG) calculates the complex amplitude scattering 
% matrix elements S1 and S2 as a function of angles ANG in radians for 
% particle size parameter X, relative index of refraction NREL. The size
% parameter X=K*R, where K is the wavenumber of incident light and R is the
% particle radius, or X=2*pi*R*NMED/WAVEL, where NMED is the refractive
% index of the medium and WAVEL is the wavelength of incident light in free
% space. The relative index of refraction NREL=NPART/NMED, where NPART is
% the index of refraction of the particle. NREL may be complex.
% 
% [S1,S2,QSCA,QEXT,QBACK] = FASTMIE(X,NREL,ANG) also calculates the
% scattering, extinction, and radar backscattering efficiencies.
%
% If ANG is scalar, then FASTMIE will generate a vector of angles for
% compatibility with the original MATLAB translation of the BHMIE Fortran
% code (Emmanuel Boss, University of Maine). In this case the vector of 
% angles is linearly spaced from 0 to pi with ANG*2-1 elements.
%
% [S1,S2,QSCA,QEXT,QBACK] = FASTMIE(X,NREL,[]) returns S1=S2=[] and skips
% calculation of all angular scattering in order to rapidly calculate
% efficiencies.
% 
% X and M must be scalar.
%
% Author: Wayne H. Slade
%         Sequoia Scientific, Inc.
%         wslade@sequoiasci.com
%
% Revisions:
% 2004/05/28 Original
% 2005/08/10 First distribution and testing (test_fastmie.m)
% 2010/05/06 Allowed empy option for nang to calculate efficiencies
% 2011/08/30 Added syntax to specify angle vector
% 2012/03/20 Documented code and streamlined calc of pi, theta

function [s1, s2, Q_scat, Q_ext, Q_back] = fastmie(x, m, nang);

    if any([length(m) length(x)]~=1)
        error('m, x must be scalar')
    end
    
	nmax = ceil(2+x+4*x.^(1/3));
	n=(1:nmax)';
	
    % Spherical Bessel functions if first (j) and second (y) and third (h, Hankel) kind. 
    % See B+H p 86 eqs 4.9-4.10 
	jx = besselj(n+0.5, x) .* sqrt(0.5*pi/x);
	jmx = besselj(n+0.5, m*x) .* sqrt(0.5*pi/(m*x));
	yx = bessely(n+0.5, x) .* sqrt(0.5*pi/x);
	hx = jx +  j*yx;
	
	j1x =  [ sin(x)/x;        jx(1:nmax-1)];
	j1mx = [ sin(m*x)/(m*x);  jmx(1:nmax-1)];
	y1x =  [ -cos(x)/x;       yx(1:nmax-1)];
	h1x = j1x + j*y1x;
    
    % Derivatives calculated using recurrence relationships, B+H p101,127
    d_jx =  x.*j1x - n.*jx;          % derivative of spherical bessel j_n(x)
	d_jmx = (m*x).*j1mx - n.*jmx;    % derivative of spherical bessel j_n(mx)
	d_hx =  x.*h1x - n.*hx;
    
    % B+H p 100, note magnetic permeability of particle = free space 
    % so mu1 = mu, which is different mu than used for calc of tau,pi
    an = ((m^2).*jmx.*d_jx - jx.*d_jmx) ./ ((m^2).*jmx.*d_hx - hx.*d_jmx);
    bn = (jmx.*d_jx - jx.*d_jmx) ./ (jmx.*d_hx - hx.*d_jmx);

    % B+H p 103
    Q_scat = (2/x.^2) * sum((2*n+1).*(abs(an).^2 + abs(bn).^2)); 
    Q_ext = (2/x.^2) * sum((2*n+1).*real(an + bn));
    % Radar backscattering efficiency, B+H p 120 section 4.6
    Q_back = (1/x.^2) * abs( sum((2*n+1).*((-1).^n).*(an-bn)) ).^2;
    
    if ~isempty(nang)
        % calculate pi and tau
        %theta = linspace(0,pi, ntheta);
        if isscalar(nang)
            theta = linspace(0,pi, fix(nang)*2-1); % calculated at (nang*2-1) angles for agreement with bhmie
        else
            theta = nang(:)';   % ensure row vector 
        end
            
        % Calculate angle-dependent functions pi and tau 
        % pi and tau will be (nmax x length(theta))
        mu = cos(theta);
        pi_n = zeros(nmax, length(theta));
        tau_n = zeros(nmax, length(theta));
        pi_n(1,:) = 1;
        pi_n(2,:) = 3*mu;
        tau_n(1,:) = mu;
        tau_n(2,:) = 2*mu.*pi_n(2,:)-3;
        for k = 3:nmax
            pi_n(k,:) = (2*k-1)/(k-1)*mu.*pi_n(k-1,:) - k/(k-1)*pi_n(k-2,:);
            tau_n(k,:) = k*mu.*pi_n(k,:) - (k+1)*pi_n(k-1,:);
        end
        
        % B+H eq 4.74 p 112
        p = repmat((2*n+1)./(n.*(n+1)), 1, length(theta)) .* pi_n;
        t = repmat((2*n+1)./(n.*(n+1)), 1, length(theta)) .* tau_n;
        AN = repmat(an, 1, length(theta));
        BN = repmat(bn, 1, length(theta));
        s1 = sum( AN.*p + BN.*t ).'; % transpose s1,s2 for dim agreement with bhmie
        s2 = sum( AN.*t + BN.*p ).';
    else
        s1 = [];
        s2 = [];
    end
    
return
    
    

