% 2D membrane clamped
% Note Courant Lambda value = sqrt(0.5)
% DTM (12/08)

% MESH DEFINITION
% f_update = c*sqrt(2)/d
% d = 10mm
% size = 1m by 1m

sr = 48000;      % sample rate of system/mesh

Nx = 10;       % number of grid points in x direction
Ny = 10;       % number of grid points in y direction
Nz = 10;       % number of grid points in z direction

Ns = 4000;      % number of samples

Inx = 5;       % Define an excitation point - can be used later
Iny = 5;       
Inz = 5;       

Outx = 5;      % Define an output point
Outy = 5;
Outz = 5;

% INITIALISE VARIABLES

% Initialise matrices used to store values at t+1, t, and t-1
p_new = zeros(Nx,Ny,Nz);       % p_new(t+1) at new time instant
p_current = p_new;      % p_new(t) at current time instant
p_previous = p_current;      % p_new(t-1) at previous time instant

% Initialise an array to capture output at a single point
out = zeros(1, Ns);

% SET EXCITATION
 
% First define a 9x9 point region in the centre of the mesh
Ne = 4;
xe = ((Nx/2)-Ne):((Nx/2)+Ne);
ye = ((Ny/2)-Ne):((Ny/2)+Ne);
ze = ((Nz/2)-Ne):((Nz/2)+Ne);
% Define a hanning window function, with an amplitude of 5.0 and a window
% size of 9 sample points - this is essenitally a smoothed impulse.
sh =  5*hanning(Ne*2+1);
% Transform this into a 2D function/shape
shape = sh*sh';
% Load the shape into the centre of the mesh at current timestep (t) 
% added two more plucking areas
p_current(xe,ye,ze) = shape;
% p_current(xe+20,ye+30) = shape;
% p_current(xe-10,ye-40) = shape;
% Let the mesh at the previous timestep, t-1 have the same shape.
p_previous = p_current;

% MESH ANIMATION

% Draw the mesh at this starting point, t=0
figure(1);
clf;
surf(p_current);
axis off;
axis equal;
shading interp;
% colormap(copper);
% V = axis;
% view([0 0])
pause;  %You need to hit a key on the keyboard to start the animation.

% set up contants
%lambda = 1 / sqrt(3);
%A = 2 - 6 * lambda^2;
%B = lambda^2;
%C = -1;
% simplified version where lamba = 1 / sqrt(3)
A = 0;
B = 1/3;
C = -1;

% MAIN LOOP
for n = (1 : Ns)
    % finite difference equation
    for k = (2 : Nx-1)
        for l = (2 : Ny-1)                              
            for m = (2 : Nz-1)   
                p_new(l, m) = A*p_current(k, l, m) ...
                            + B*(p_current(k+1, l, m) + p_current(k-1, l, m) ...
                               + p_current(k, l+1, m) + p_current(k, l-1, m) ...
                               + p_current(k, l, m+1) + p_current(k, l, m-1)) ...
                            + C*p_previous(k, l, m);
            end
        end
    end
    
    
    % Get Output Value
    
    out(n) = p_new(Outx, Outy, Outz);
   
    % Update mesh history
    
    p_previous = p_current;
    p_current = p_new;
    
    % update plot
    surf(p_new);
    axis off;
    axis equal;
%     axis(V);
    shading interp;
%     view([0 0])
    pause(0.001);
end

% Plot time domain response of output point
figure(2);
clf;
plot(out);
xlabel('time (samples)');
ylabel('Amplitude');

fftSize = 8192;
f = (0:fftSize-1)*(sr/fftSize);

% Plot frequency domain response of output point
figure(3);
clf;
semilogx(f,20*log10(abs(fft(out,fftSize))/max(abs(fft(out,fftSize)))));
axis([0 8000 -40 0]);
xlabel('frequency (Hz)');
ylabel('magnitude response (dB)');

h = fir1(20,0.25); % 20th order low-pass filter with cutoff at 
                   % 0.25 times the sample rate normalised to half-Nyquist
lowpassout = filter(h,1,out); % filtered version of ‘out’ using ‘h’
sound(lowpassout, sr, 16); % plays a vector as a sound
