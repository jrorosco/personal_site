%PLATINUM_MODEL Compute electrical, optical, and radiative properties of platinum.
%
% Given the relevant temperature, 'T', wavelength, 'l', and/or incident (or
% emission) angle, 'theta_i', returns the property indicated by the type
% parameter, 'type'.
%
%   Inputs:
%
%       l       : wavelength [um], double
%       T       : temperature [K], double
%       theta_i : incident (or emission) angle [rad], double
%       type    : property type indicator, string
%                   'su' : susceptibility
%                   'pe' : permittivity
%                   'nk' : refractive index
%                   'em' : emissivity
%                   're' : reflectivity
%                   'dc' : dc resistivity
%
%   Outputs:
%
%       ret_val : depends on type selected by input
%                   'su' : susceptibility, complex double
%                   'pe' : permittivity, complex double
%                   'nk' : refractive index, complex double
%                   'em' : emissivity, real double
%                   're' : reflectivity, real double
%                   'dc' : dc resistivity [ohm*um], real double
%
%   Notes:
%
%       (1) Relies on Faddeeva package implemented by Steven G. Johnson,
%           available at http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package
%
%   Examples of usage:
%
%       normal_emissivity = platinum_model(8,1000,0,'em');
%       dc_resistivity = platinum_model([],1000,[],'dc');
%       refractive_index = platinum_model(linspace(1.5,16,100),1000,45,'nk');

% J. Orosco and C. F. M. Coimbra
% Coimbra Research Group
% Department of Mechanical
% and Aerospace Engineering
% University of California, San Diego

% Copyright 2019 Jeremy Orosco
% This work is licensed under the
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0
% International License. To view a copy of this license,
% visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function ret_val = platinum_model(l,T,theta_i,type)

    % physical constants
    c0 = 299792458;             % vac. light speed          [m/s]
    e0 = 8.85418782*10^(-12);   % vac. permittivity         [s^4*A^2*m^-3*kg^-1]
    T_ref = 294;                % reference temperature     [K]
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% resistivity and conductivity %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % high-temp. resistivity parms
    rho_h = 0.1137;             % power-law multiplier      [ohm*um]
    k_h = 0.9197;               % power law                 [null]
    
    % low-temp. resistivity parms
    rho_0 = 0.0928;             % residual resistivity      [ohm*um]
    rho_l = rho_h-rho_0;        % power-law multiplier      [ohm*um]
    k_l = k_h*(rho_h/rho_l);    % power law                 [null]
    
    % low-temp. resistivity
    rho_low =@(T) rho_l*(T/T_ref).^k_l+rho_0;
    
    % high-temp. resistivity
    rho_high =@(T) rho_h*(T/T_ref).^k_h;

    % conditional resistivity
    rho = (T/T_ref>=1).*rho_high(T) + (T/T_ref<1).*rho_low(T);
    
    % ratiometric parms
    kappa = 0.4590;             % auxiliary parameter       [uohm/um]
    sigma_a = 0.2530;           % anomalous conductivity    [ohm*um]
    
    % conductivity
    sigma = 1./rho;
    lambda = sigma/kappa;
    alpha = sigma/sigma_a;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%% material parameters / functions %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % interband model parms
    x1 = sigma;             % conductivity                  [ohm*um]
    x2 = lambda;            % relaxation wavelength         [um]
    x3 = 18.2500;           % disorder wavelength           [um]
    x4 = 0.0982;            % memory decay strength         [null]
    x5 = alpha;             % time constant shift           [null]
    x6 = 0.5800;            % intra. oscillator strength    [null]
    
    % intraband model parms
    x7 = 0.1293;            % plasma wavelength ref.        [um]
    x8 = 0.0026;            % nondim. coeff. of therm. exp. [null]
    x9 = 2.4533;            % damping wavelength ref.       [um]
    x10 = 0.4012;           % damping wavelength shift      [null]
    x11 = 1.5457;           % 0K resonant wavelength        [um]
    x12 = 0.0365;           % Varshni shift                 [null]
    x13 = 0.8019;           % Varshni cutoff                [null]
    x14 = 8.0450;           % gaussian disorder wavelength  [um]
    x15 = 0.1887;           % inter. oscillator strength    [null]
    
    % thermal expansion factor for electron density
    k_therm = 1./(1+x8*(T/T_ref-1));
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% intraband component %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % drude susceptbility
    chi_d = -(x1*x6.*k_therm.*l.^2)./((2.*pi.*c0.*e0).*(x2+1j*l));
    
    % anomalous susceptibility
    chi_a = -((-1j).^x4).*(l.^2.*(x1./x5)*x6.*k_therm.*(x3./l).^x4)./ ...
            ((2.*pi.*c0.*e0).*(x2./x5+1j.*l));
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% interband component %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % small parm
    delta = sqrt(eps);
    
    % T-dep parms
    lp = x7./sqrt(k_therm);                                 % plasma wavelength
    ld = (x9.^-1+x10*(sqrt(T/T_ref)-1)).^(-1);              % damping wavelength
    
    % wavelength-dependent Varshni equation
    lr = (x11^-1-x12*(T/T_ref).^2./((T/T_ref)+x13)).^-1;    % resonant wavelength

    % alpha parms
    alpha1 = ((2*l).^(-1/2)).*((l.^-2+ld.^-2).^(1/2)+l.^-1).^(1/2);
    alpha2 = ((2*l).^(-1/2)).*((l.^-2+ld.^-2).^(1/2)-l.^-1).^(1/2)+delta;
    alpha = alpha1 + 1j*alpha2;

    % auxiliary parms
    zp = (alpha - lr.^-1)./(sqrt(2)*x14.^-1);
    zm = (-alpha - lr.^-1)./(sqrt(2)*x14.^-1);

    % auxiliary exp-log product funs
    felp = exp(-zp.^2).*(log(zp)+log(-1./zp)-1j*pi);
    felm = exp(-zm.^2).*(log(zm)+log(-1./zm)-1j*pi);

    % auxiliary model funs
    fwp = 1j*pi*Faddeeva_w(zp) + felp;
    fwm = 1j*pi*Faddeeva_w(zm) + felm;

    % normalization factor
    chi0 = -4*sqrt(pi)*Faddeeva_Dawson(-x14./(sqrt(2)*lr));

    % dc amplitude
    A = (lr.^2*x15)./(lp.^2);

    % shape function
    S = (fwp+fwm)./abs(chi0);

    % oscillator
    chi_b = A.*S;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% optical properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % total susceptibility
    chi = chi_d + chi_a + chi_b;
    
    % permittivity
    K = 1 + chi;
    
    % refractive index
    N = sqrt(K);
    
    % indices
    n = real(N);
    k = imag(N);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% radiative properties %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % auxiliary fresnel equations
    nt2 = n.^2-k.^2-sin(theta_i).^2;
    pf = sqrt((1/2)*(sqrt(nt2.^2+(2*n.*k).^2)+nt2));
    qf = sqrt((1/2)*(sqrt(nt2.^2+(2*n.*k).^2)-nt2));
    
    % fresnel s-polarized reflection coefficient
    rsn = cos(theta_i)-pf+1j*qf;
    rsd = cos(theta_i)+pf-1j*qf;
    rs = rsn./rsd;
    
    % fresnel p-polarized reflection coefficient
    rpn = pf-sin(theta_i).*tan(theta_i)-1j*qf;
    rpd = pf+sin(theta_i).*tan(theta_i)-1j*qf;
    rp = rs.*(rpn./rpd);
    
    % polarized reflectivity
    Rs = abs(rs).^2;
    Rp = abs(rp).^2;
    
    % unpolarized reflectivity
    re = (1/2)*(Rs+Rp);
    
    % unpolarized emissivity
    em = 1-re;
    
    % set return value
    switch type
        case 'su'
            ret_val = chi;
        case 'pe'
            ret_val = K;
        case 'nk'
            ret_val = N;
        case 'em'
            ret_val = em;
        case 're'
            ret_val = re;
        case 'dc'
            ret_val = rho;
    end

end