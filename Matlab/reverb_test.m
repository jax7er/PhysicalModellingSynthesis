function reverb_test

%   Copyright 2011 The MathWorks, Inc.

    % Load sample audio data (s.y=audio data, s.fs=SampleRate)
    s = load('speech_dft.mat');

    % Apply reverberbation effect on audio data
    tic;
    runMexVersion = exist('do_reverb_mex','file') == 3;
    % Is the MEX version available?
    if runMexVersion
        iter = 100;
        for i = 1:iter
            y = do_reverb_mex(s.y,s.fs);
        end
    else
        % Default on executing function in standard MATLAB
        y = do_reverb(s.y,s.fs);
    end
    dt = toc;
    if runMexVersion
        dt = dt/iter;
    end
    
    if audiodevinfo(0, s.fs, 16, 1) ~= -1 % Is audio output available?
        % Play result on audio output
        p = audioplayer(y, s.fs);
        playblocking(p);
    end
    
    fprintf('Running time = %d milliseconds\n',round(dt*1000));
end
