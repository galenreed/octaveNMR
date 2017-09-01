clear all;
close all;


numParticles = 2;
stateSpaceMatrixSize = 2^numParticles;
J = 5; % [Hz]
sx = [0 1;
      1 0]/2;
sy = [0 -1i;
      1i 0]/2;
sz = [1  0;
      0 -1]/2;      
si = [1  0;
      0  1]/2;

% define a bunch of density matrix elements
I1x = kron(sx,si) * 2;
I2x = kron(si,sx) * 2;
I1y = kron(sy,si) * 2;
I2y = kron(si,sy) * 2;
I1z = kron(sz,si) * 2;
I2z = kron(si,sz) * 2;

% Bilinear operators; these are actually 2I1zI2x etc ... 
I1zI2z = kron(sz,sz) * 2; 
I1zI2x =  kron(sz,sx) * 2;
I1zI2y =  kron(sz,sy) * 2;


% pulse sequence. all vectors should be the same length 
% for the RF, put angles and corresponding delays of 0
% for delay periods, set RF to 0     
tau = 1/(2*J);

% refocused inept
channel1RF = [90 0 180 0 90*1i] * pi/180;
channel2RF = [0  0 180 0 90] * pi/180;
delays = [0 tau/2 0 tau/2 0];

% basic inept
%channel1RF = [90 0 90*1i] * pi/180;
%channel2RF = [0  0 90] * pi/180;
%delays =     [0 tau 0];



driftHamiltonian = pi * J * I1zI2z; % J coupling and chemical shift, [Hz]
numPulseSequenceElements = length(delays);
rho = I1z;% initial state
densityMatrix = zeros([stateSpaceMatrixSize stateSpaceMatrixSize numPulseSequenceElements]);

for ii = 1:numPulseSequenceElements
  currentHamiltonian = [];
  if(delays(ii) == 0) % rf hamiltonian
    RF1x = real(channel1RF(ii));
    RF1y = imag(channel1RF(ii));
    RF2x = real(channel2RF(ii));
    RF2y = imag(channel2RF(ii));
    currentHamiltonian = RF1x*I1x +  RF1y*I1y + RF2x*I2x + RF2y*I2y;
  else % drift hamiltonian
    currentHamiltonian = driftHamiltonian * delays(ii);
  end
  Ht = currentHamiltonian;
  rho = expm(-1i * Ht) * rho * expm(+1i * Ht);
  densityMatrix(:,:,ii) = rho;
end



% frobenius inner product
function innerProduct = fip(A,B)
  innerProduct = trace(conj(A') * B);
endfunction


% track spin 2 density matrix elements over pulse sequence
spin2X = zeros([1 numPulseSequenceElements]);
spin2Y = zeros([1 numPulseSequenceElements]);
spin1Zspin2X = zeros([1 numPulseSequenceElements]);
spin1Zspin2Y = zeros([1 numPulseSequenceElements]);

for ii = 1:numPulseSequenceElements
  rho = squeeze(densityMatrix(:,:,ii));
  spin2X(ii) = fip(rho, I2x);
  spin2Y(ii) = fip(rho, I2y);
  spin1Zspin2X(ii) = fip(rho, I1zI2x);
  spin1Zspin2Y(ii) = fip(rho, I1zI2y);
end

t = 1:numPulseSequenceElements;
plot(t,spin2X,...
     t,spin2Y,...
     t,spin1Zspin2X,...
     t,spin1Zspin2Y)
legend('I2x','I2y','2I1zI2x','2I1zI2y');








     