% John P. Sutter
% 28 June 2011
% This program calculates the final wave vector and final electric field
% amplitudes for a given incident wave vector with initial electric field
% amplitudes. Dynamical diffraction theory for a perfect crystal is 
% applied. Only the 2-beam Bragg case with diffraction plane normal to
% crystal surface is considered. Equations are from William H. Zachariasen,
% "Theory of X-Ray Diffraction in Crystals" Chapter III.
%
% CONSTANTS
% For converting structure factors into electric susceptibility Fourier
% components
epmcA = 8.969785162E-06;
% Silicon lattice parameters
% 300 K
aSiART = 5.43108124;
% 120 K: minimum value
aSiALT = 5.42965775;
% Germanium lattice parameters (for 70Ge):
% 300 K
aGeART = 5.6579830;
% 75 K
aGeALT = 5.6525930;
% Diamond lattice parameters
% 298 K
aDiART = 3.56712;
% 77 K
aDiALT = 3.56681;
% Sapphire lattice parameters (hexagonal axes)
% 300 K
aSaART = 4.7593240;
cSaART = 12.9919261;
% 80 K
aSaALT = 4.7563257;
cSaALT = 12.9820703;
% Quartz lattice parameters (hexagonal axes)
% 296 K
aQuART = 4.9141;
cQuART = 5.4060;
% 78 K
aQuALT = 4.9030;
cQuALT = 5.3999;
% This is different for every crystal.
% These crystal parameters are tabulated. A good source of information 
% on crystal structures (lattice parameters and atomic positions in the
% unit cell) is the set of volumes by R. Wyckoff.
% END OF CONSTANTS
%
% USER INPUTS
% Input #0a: crystal material
xtltyp = 1;
% Input #0b: choice of room-temperature or low-temperature lattice
% parameters: rmtmp = 0 for low temperature, 1 for high temperature.
rmtmp = 1;
% Calculations from input #0a,b: volume of crystal lattice unit cell
if isequal(xtltyp,1)
    % silicon
    if isequal(rmtmp,0)
        aSiA = aSiALT;
    else
        aSiA = aSiART;
    end
    volA3 = aSiA^3;
elseif isequal(xtltyp,2)
    % germanium
    if isequal(rmtmp,0)
        aGeA = aGeALT;
    else
        aGeA = aGeART;
    end
    volA3 = aGeA^3;
elseif isequal(xtltyp,3)
    % diamond
    if isequal(rmtmp,0)
        aDiA = aDiALT;
    else
        aDiA = aDiART;
    end    
    volA3 = aDiA^3;
elseif isequal(xtltyp,4)
    % sapphire
    if isequal(rmtmp,0)
        aSaA = aSaALT;
        cSaA = cSaALT;
    else
        aSaA = aSaART;
        cSaA = cSaART;
    end    
    volA3 = sqrt(3.)*aSaA^2*cSaA/2.;
elseif isequal(xtltyp,5)
    % quartz
    if isequal(rmtmp,0)
        aQuA = aQuALT;
        cQuA = cQuALT;
    else
        aQuA = aQuART;
        cQuA = cQuART;
    end        
    volA3 = sqrt(3.)*aQuA^2*cQuA/2.;
end
% Input #1: photon energy in eV
EceV = 11183.9195;
% Calculations from input #1
alamA = 12398.4193009/EceV;
k0Ai = 1./alamA;
f0c = 113.3970+1.3727i;
% structure factor of O beam. Structure factors are not given directly in
% Wyckoff. However, Wyckoff gives the atomic positions inside the unit
% cell. With this information, one can calculate the structure factors for
% any reflection using a simple formula.
psi0c = -epmcA*alamA^2*f0c/volA3
% 0th Fourier component of the crystal's electric susceptibility.
% Input #2: Miller indices of Bragg reflection. See Kittel for definitions.
hMilND = 4;
kMilND = 4;
lMilND = 4;
% Calculations from input #2: interplanar spacing (note dependence on
% input #0!)
if isequal(xtltyp,1)
    % silicon
    dA = aSiA/sqrt(hMilND^2 + kMilND^2 + lMilND^2);
elseif isequal(xtltyp,2)
    % germanium
    dA = aGeA/sqrt(hMilND^2 + kMilND^2 + lMilND^2);
elseif isequal(xtltyp,3)
    % diamond
    dA = aDiA/sqrt(hMilND^2 + kMilND^2 + lMilND^2);
elseif isequal(xtltyp,4)
    % sapphire
    dA = 1./sqrt( 4.*(hMilND^2+kMilND^2+hMilND*kMilND)/(3.*aSaA^2) + lMilND^2/cSaA^2 );
elseif isequal(xtltyp,5)
    % quartz
    dA = 1./sqrt( 4.*(hMilND^2+kMilND^2+hMilND*kMilND)/(3.*aQuA^2) + lMilND^2/cQuA^2 );
end
fhc = -34.4013-1.1447i; % structure factor for diffracted beam
psihc = -epmcA*alamA^2*fhc/volA3
fmhc = -59.6881-0.7831i; % structure factor for diffracted->incident bounce
psimhc = -epmcA*alamA^2*fmhc/volA3
% Input #3: In-plane asymmetry angle:
% This angle is POSITIVE for grazing incidence and NEGATIVE for grazing
% exit!
alphdg = 0.;
% Calculations from input #3
alphrd = alphdg*pi/180.;
HXAi = [ 0.; cos(alphrd)/dA; -sin(alphrd)/dA ]
% Input #4: Input mode type. Two modes are available:
% AUTO: User inputs the deviation from the Bragg angle.
% ABSOLUTE: User inputs the absolute angle relative to the diffracting
% planes. This is not equal to the pitch rotation angle of the crystal if
% the asymmetry angle is nonzero.
ptMode = 0;
if isequal(ptMode,0)
    % AUTO mode: Bragg angle is calculated automatically.
    % Input #4a: energy EadjeV at which crystal's Bragg condition should be
    % fulfilled when dthur (input #4b) is zero. Note: EadjeV need not be
    % the central energy of the wavefront!
    EadjeV = 11183.9195;
    % Input #4b: deviation dthur from kinematic Bragg angle at E = EadjeV
    dthur = 0.;
    % Calculations from input #4a
    aladjA = 12398.4193009/EadjeV;
    thBrd = asin(aladjA/(2.*dA));
    thrd = thBrd + dthur*1.e-06;
elseif isequal(ptMode,1)
    % ABSOLUTE mode: absolute value of incidence angle is entered by user
    % Input #4c: absolute incidence angle from diffracting planes.
    thdg = 42.;
    % Calculations from input #4c
    thrd = thdg*pi/180.;
end
% Calculations from input #4: thptrd = the pitch rotation angle of the
% crystal in radians. This is the angle of rotation of the crystal
% coordinate system relative to the lab frame.
thptrd = thrd - alphrd
% Input #5: Roll angle
chidg = 0.;
% Calculations from input #6
chird = chidg*pi/180.;
% Input #6: Yaw angle
phidg = 0.;
% Calculations from input #7
phird = phidg*pi/180.;
% Input #7: Thickness in microns
thicum = 500.;
% Input #8: Whether to calculate the transmitted beam as well as the
% diffracted beam. itrans = 0 NO, itrans = 1 YES.
itrans = 1;
% END OF USER INPUTS
%
% TRANSFORMATION MATRIX: 
% RXtLab: 3x3 orthogonal matrix that converts 3x1 vector components in
% crystal coordinates to components in lab coordinates
% RLabXt: transpose of RXtLab: converts 3x1 vector components in lab
% coordinates to components in crystal coordinates
RXtLab(1,1) =  cos(chird)*cos(phird);
RXtLab(1,2) = -sin(chird);
RXtLab(1,3) =  cos(chird)*sin(phird);
RXtLab(2,1) =  cos(thptrd)*sin(chird)*cos(phird)-sin(thptrd)*sin(phird);
RXtLab(2,2) =  cos(thptrd)*cos(chird);
RXtLab(2,3) =  cos(thptrd)*sin(chird)*sin(phird)+sin(thptrd)*cos(phird);
RXtLab(3,1) = -sin(thptrd)*sin(chird)*cos(phird)-cos(thptrd)*sin(phird);
RXtLab(3,2) = -sin(thptrd)*cos(chird);
RXtLab(3,3) = -sin(thptrd)*sin(chird)*sin(phird)+cos(thptrd)*cos(phird);
RLabXt = RXtLab';
%
% DYNAMICAL DIFFRACTION CALCULATIONS
% This is a loop over various values of (kx,ky).
nk = 1001;
empty = [];
save TestK.out empty -ascii;
save TestRef.out empty -ascii;
save TestAmp.out empty -ascii;
save TestTrns.out empty -ascii;
% Polarization vectors e1, e2 in lab frame:
% z along beam
% y upward
% x = y x z
% xz-plane is horizontal; yz-plane is vertical
% NOTE: We neglect dependence of polarization vectors on (kx,ky). If the
% angular range of the wave vectors is under 1 mrad, this should be OK.
e1 = [ 1.; 0.; 0. ];
e2 = [ 0.; 1.; 0. ];
% Conversion of polarization vector components to crystal frame:
e1X = RLabXt*e1;
e2X = RLabXt*e2;
% Calculation of new polarization vectors in the crystal frame, with the
% sigma vector parallel to the diffracting atomic planes of the crystal. In
% this simple case, the sigma vector will also be parallel to the crystal
% surface.
% This is an especially convenient choice because the sigma and pi
% components of the wavefields inside the crystal can be treated
% separately.
% Again, we neglect the dependence on (kx,ky).
kc0XAi = RLabXt*[0.;0.;k0Ai];
ucn0X = kc0XAi/k0Ai;
sg0X = cross(HXAi,ucn0X);
sg0X = sg0X/norm(sg0X);
pi0X = cross(ucn0X,sg0X);
% Calculate transformation matrix from (e1X,e2X) to (sg0X,pi0X)
PolTrn = [ dot(e1X,sg0X), dot(e2X,sg0X); dot(e1X,pi0X), dot(e2X,pi0X) ];
% Calculate polarization vectors in the crystal frame for the diffracted
% beam, ignoring the dependence on (kx,ky):
HtXAi = [HXAi(1); 0.; HXAi(3)];
kcHtXA = [ kc0XAi(1); 0.; kc0XAi(3) ] + HtXAi;
kcHXAi = kcHtXA + [0.; sqrt(k0Ai^2-dot(kcHtXA,kcHtXA)); 0.];
ucnHX = kcHXAi/norm(kcHXAi);
sgHX = sg0X;
piHX = cross(ucnHX,sgHX);
% Transform the diffracted beam polarization vectors back into the lab
% frame:
sgH = RXtLab*sgHX;
piH = RXtLab*piHX;
% Begin loop over wave vectors:
for ik = 1:nk
% Hard-coded inputs: wave vector components
   kxAi =  0.00; %inverse Angstroms
   kyAi = -0.00010+0.00020*(ik-1)/(nk-1); %inverse Angstroms
% Hard-coded inputs: electric field amplitudes in k-space
   EIn12C = [ 0.+0.i ; 1.+0.i ];
% Ein12C(1) along e1
% Ein12C(2) along e2
% k0wvAi = incident beam wave vector (kx,ky,kz).
   kzAi = sqrt(k0Ai^2-kxAi^2-kyAi^2);
   k0wAi = [kxAi;kyAi;kzAi];
   u0 = k0wAi/norm(k0wAi);
% Conversion of wave vector components to crystal frame:
   k0wXAi = RLabXt*k0wAi;
   u0X = k0wXAi/norm(k0wXAi);
% Convert components of incident beam polarization from (e1X,e2X) to
% (sg0X,pi0X).
   EInSPC = PolTrn*EIn12C;
% Calculate the wave vector for the diffracted beam.
% Note that index of refraction corrections are included.
   k0tXAi = [k0wXAi(1); 0.; k0wXAi(3)];
   kHtXAi = k0tXAi + HtXAi;
   kmHXAi = kHtXAi + [0.; sqrt(k0Ai^2-dot(kHtXAi,kHtXAi)); 0.];
% Calculate direction cosine gamma0, reflection asymmetry parameter bee, 
% deviation parameter Adev and normalized deviation parameter zeeC (as 
% defined by Zachariasen gamma0, b, alpha and z).
   gamma0 = -u0X(2);
   bee = 1./(1.+HXAi(2)/k0wXAi(2));
   Adev = ( 2.*dot(k0wXAi,HXAi) + 1./dA^2 )/k0Ai^2;
   zeeC = 0.5*( (1.-bee)*psi0c + bee*Adev );
% Calculate the complex reflectivity DHsgC for sigma polarization.
   queC = bee*psihc*psimhc;
   sqrqzC = sqrt(queC+zeeC^2);
   del1C = 0.5*( psi0c - zeeC + sqrqzC );
   del2C = 0.5*( psi0c - zeeC - sqrqzC );
   x1C = (-zeeC+sqrqzC)/psimhc;
   x2C = (-zeeC-sqrqzC)/psimhc;
   Cph1C = exp(-2.*pi*1.i*k0Ai*del1C*(thicum*1.E+04)/gamma0);
   Cph2C = exp(-2.*pi*1.i*k0Ai*del2C*(thicum*1.E+04)/gamma0);
   if (isinf(Cph1C))
       DHsgC = x2C;
   elseif (isinf(Cph2C))
       DHsgC = x1C;
   else
       DHsgC = x1C*x2C*(Cph2C-Cph1C)/(Cph2C*x2C-Cph1C*x1C);
   end
   if isequal(itrans,1)
       % calculate the complex reflectivity D0trsC of the transmitted beam.
       if (isinf(Cph1C))
           D0trsC = -Cph2C*(x2C-x1C)/x1C;
       elseif (isinf(Cph2C))
           D0trsC = +Cph1C*(x2C-x1C)/x2C;
       else
           D0trsC = Cph1C*Cph2C*(x2C-x1C)/(Cph2C*x2C-Cph1C*x1C);
       end
   end
% Calculate the complex reflectivity DHpiC for pi polarization.
   queC = bee*psihc*psimhc*(cos(2.*thrd))^2;
   sqrqzC = sqrt(queC+zeeC^2);
   del1C = 0.5*( psi0c - zeeC + sqrqzC );
   del2C = 0.5*( psi0c - zeeC - sqrqzC );
   x1C = (-zeeC+sqrqzC)/(psimhc*cos(2.*thrd));
   x2C = (-zeeC-sqrqzC)/(psimhc*cos(2.*thrd));
   Cph1C = exp(-2.*pi*1.i*k0Ai*del1C*(thicum*1.E+04)/gamma0);
   Cph2C = exp(-2.*pi*1.i*k0Ai*del2C*(thicum*1.E+04)/gamma0);
   if (isinf(Cph1C))
       DHpiC = x2C;
   elseif (isinf(Cph2C))
       DHpiC = x1C;
   else
       DHpiC = x1C*x2C*(Cph2C-Cph1C)/(Cph2C*x2C-Cph1C*x1C);
   end 
   if isequal(itrans,1)
       % calculate the complex reflective D0trpC of the transmitted beam.
       if (isinf(Cph1C))
           D0trpC = -Cph2C*(x2C-x1C)/x1C;
       elseif (isinf(Cph2C))
           D0trpC = +Cph1C*(x2C-x1C)/x2C;
       else
           D0trpC = Cph1C*Cph2C*(x2C-x1C)/(Cph2C*x2C-Cph1C*x1C);
       end
   end
% Calculate the diffracted amplitudes:
   EHSPC = [DHsgC,0.;0.,DHpiC]*EInSPC;
   if isequal(itrans,1)
       E0tSPC = [D0trsC,0.;0.,D0trpC]*EInSPC;
   end
% Convert the diffracted beam's wave vector from the crystal frame back 
% to the lab frame:
   kmHAi = RXtLab*kmHXAi;
% for test output only
   dangur = -asin(kyAi/k0Ai)*1.E+06;
   uref = [ 0. ; sin(2.*thBrd) ; cos(2.*thBrd) ];
   kangur = acos(dot(uref,kmHAi)/k0Ai)*1.E+06;
   if (abs(imag(kangur)) > 0.)
       kangur = 0.;
   end
   ukcrAi = cross(uref,kmHAi);
   if (ukcrAi(1) > 0.)
       kangur = -kangur;
   end
   kmat = [ dangur, norm(kmHAi)/k0Ai, kangur ];
   save TestK.out kmat -append -ascii
   refmat = [dangur,abs(DHsgC)^2,angle(DHsgC),abs(DHpiC)^2,angle(DHpiC)];
   save TestRef.out refmat -append -ascii
   ampmat = [dangur,abs(EHSPC(1))^2,angle(EHSPC(1)),abs(EHSPC(2))^2,angle(EHSPC(2))];
   save TestAmp.out ampmat -append -ascii
   if isequal(itrans,1)
       trmat = [dangur,abs(E0tSPC(1))^2,angle(E0tSPC(1)),abs(E0tSPC(2))^2,angle(E0tSPC(2))];
       save TestTrns.out trmat -append -ascii
   end
end