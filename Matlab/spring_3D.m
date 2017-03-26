% 3D solid clamped
% Note Courant Lambda value = sqrt(1/3)

% MESH DEFINITION
% f_update = c*sqrt(2)/d

sr = 44100;      % sample rate of system/mesh
d = 330*sqrt(2)/sr;

Nx = 100;       % number of grid points in x direction
Ny = 100;       % number of grid points in y direction
Nz = 100;       % number of grid points in z direction

fprintf('size_x = %f m\n', d*Nx);
fprintf('size_y = %f m\n', d*Ny);
fprintf('size_z = %f m\n', d*Nz);

seconds = 2;   % 10s takes roughly 8 hours to calculate with 100x100x100 and 44.1kHz
Ns = sr*seconds;      % number of samples

Inx = Nx/2;       % Define an excitation point - can be used later
Iny = Ny/2;       
Inz = Nz/2;       

Outx = Nx/2;      % Define an output point
Outy = Ny/2;
Outz = Nz/2;

% INITIALISE VARIABLES

% SET EXCITATION
 
% First define an Se by Se point region in the centre of the mesh
Ne = 2;
Se = Ne*2+1;
% Define a hanning window function, with an amplitude of 5.0 and a window
% size of Se sample points - this is essenitally a smoothed impulse.
han = hanning(Se);
han_2d = han*han';
han_3d = zeros(Se, Se, Se);
for z_count = (1 : Se)
    han_3d(:,:,z_count) = han_2d * han(z_count);
end
han_3d = han_3d * 3;

% Initialise matrices used to store values at t+1, t, and t-1

p_current = [zeros(Nx,Ny,Nz) zeros(Nx,Ny,Nz) zeros(Nx,Ny,Nz)];       % p_new(t+1) at new time instant
% p_new = [zeros(Nx,Ny,Nz) zeros(Nx,Ny,Nz) zeros(Nx,Ny,Nz)];       % p_new(t+1) at new time instant
% p_current = p_new;      % p_new(t) at current time instant
% p_previous = p_current;      % p_new(t-1) at previous time instant

px_new = [zeros(Nx,Ny,Nz) zeros(Nx,Ny,Nz) zeros(Nx,Ny,Nz)];       % p_new(t+1) at new time instant
px_current = px_new;      % p_new(t) at current time instant
px_previous = px_current;      % p_new(t-1) at previous time instant

py_new = [zeros(Nx,Ny,Nz) zeros(Nx,Ny,Nz) zeros(Nx,Ny,Nz)];       % p_new(t+1) at new time instant
py_current = py_new;      % p_new(t) at current time instant
py_previous = py_current;      % p_new(t-1) at previous time instant

pz_new = [zeros(Nx,Ny,Nz) zeros(Nx,Ny,Nz) zeros(Nx,Ny,Nz)];       % p_new(t+1) at new time instant
pz_current = pz_new;      % p_new(t) at current time instant
pz_previous = pz_current;      % p_new(t-1) at previous time instant

% Initialise an array to capture output at a single point
out = zeros(1, Ns);

xe_range = (Inx-Ne):(Inx+Ne);
ye_range = (Iny-Ne):(Iny+Ne);
ze_range = (Inz-Ne):(Inz+Ne);
% Load the shape into the centre of the mesh at current timestep (t) 
% added two more plucking areas
px_current(xe_range, ye_range, ze_range) = han_3d;
% Let the mesh at the previous timestep, t-1 have the same shape.
px_previous = px_current;
py_current(xe_range, ye_range, ze_range) = han_3d;
% Let the mesh at the previous timestep, t-1 have the same shape.
py_previous = py_current;
pz_current(xe_range, ye_range, ze_range) = han_3d;
% Let the mesh at the previous timestep, t-1 have the same shape.
pz_previous = pz_current;

% set up contants
% scale = 0.99999;
scale = 1;

courant = 0;
Bcourant = scale / 3;
Ccourant = -scale;

if courant == 1
    lambda_squared = 1 / sqrt(3);
else
    lambda_squared = 1 / 5;
end
A = scale * (2 - 6 * lambda_squared);
B = scale * (lambda_squared);
C = scale * (-1);

graph = 1;
figure(1);
xyones = ones(Nx, Ny);
filename = sprintf('FD-3D Nx_%d Ny_%d Nz_%d d_%f Ne_%d lambda_%f scale_%f sr_%d len_%d ts_%s', Nx, Ny, Nz, d, Ne, sqrt(lambda_squared), scale, sr, seconds, replace(replace(datestr(now), ':', '-'), ' ', '_'));
dispstat('','init');

spring_const = 0.1;
damp_fact = 0.1;
spacing = 0.1;
mass = 0.001;
% const_scale = -(spring_const+(damp_fact/spacing))/(2*mass);
const_scale = 0.1;

% MAIN LOOP
for n = (1 : Ns)
    % finite difference equation
    for k = (2 : Nx-1)
        for l = (2 : Ny-1)                              
            for m = (2 : Nz-1)
                px_new(k, l, m) = (2 - 2 * const_scale) * px_current(k, l, m) ...
                               + const_scale * (px_current(k+1, l, m) + px_current(k-1, l, m)) ...
                               - px_previous(k, l, m);
                py_new(k, l, m) = (2 - 2 * const_scale) * py_current(k, l, m) ...
                               + const_scale * (py_current(k, l+1, m) + py_current(k, l-1, m)) ...
                               - py_previous(k, l, m);
                pz_new(k, l, m) = (2 - 2 * const_scale) * pz_current(k, l, m) ...
                               + const_scale * (pz_current(k, l, m+1) + pz_current(k, l, m-1)) ...
                               - pz_previous(k, l, m);
                
                p_current(k, l, m) = sqrt(px_new(k, l, m)^2 + py_new(k, l, m)^2 + pz_new(k, l, m)^2);
            end % m
        end % l
    end % k
        
    % Get Output Value
    
    out(n) = p_current(Outx, Outy, Outz);
   
    % Update mesh history
    
    px_previous = px_current;
    px_current = px_new;
    py_previous = py_current;
    py_current = py_new;
    pz_previous = pz_current;
    pz_current = pz_new;
    
    % update plot
    if not(mod(n/Ns, 0.01))
        dispstat(sprintf('%d%%',100*n/Ns));
    end
        
    if graph == 1
        clf;
%         for b = (2 : Nz-1)
%             mesh(p_current(:,:,b) + b, 'facecolor', 'none');
%             hold on
%         end
        surface(p_current(:,:,Outz) + Outz);
        axis equal;
        axis off;
        axis([1 Nx 1 Ny 1 Nz]);
        shading interp;
%         view([0 0]);
        pause(0.01);
    end
end

% Plot time domain response of output point
figure(2);
clf;
plot(out);
xlabel('time (samples)');
ylabel('Amplitude');
% time_filename = sprintf('output/time_%s.fig', filename);
% saveas(figure(2), time_filename);

fftSize = 8192;
f = (0:fftSize-1)*(sr/fftSize);

% Plot frequency domain response of output point
figure(3);
clf;
semilogx(f,20*log10(abs(fft(out,fftSize))/max(abs(fft(out,fftSize)))));
axis([0 8000 -40 0]);
xlabel('frequency (Hz)');
ylabel('magnitude response (dB)');
% frequency_filename = sprintf('output/frequency_%s.fig', filename);
% saveas(figure(2), frequency_filename);

h = fir1(20,0.25); % 20th order low-pass filter with cutoff at 
                   % 0.25 times the sample rate normalised to half-Nyquist
lowpassout = filter(h,1,out) ./ max(abs(out)); % filtered version of ‘out’ using ‘h’
sound(lowpassout, sr, 16); % plays a vector as a sound

% audio_filename = sprintf('output/audio_%s.wav', filename);
% audiowrite(audio_filename, lowpassout, sr);