% Convolutional reverb - Offline
% T. Hueber - CNRS/GIPSA-lab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

% read the sample waveform
filename='./MG_list1_sent381.wav';
[x,Fs] = audioread(filename);

buffer_size = 512;
% read the impulse response waveform
filename='./IMreverbs1/Five Columns.wav';
[imp,Fsimp] = audioread(filename);

% Keep only one channel 
imp_left = imp(:,1);

M = length(imp_left);
L = length(x);

step_fft =1;
overlap_add = 1;

if (step_fft == 0)
  % convolution in the temporal domain (slower)
  if(overlap_add == 0)
      x_conv = zeros(1,L+M-1);
      for n=1:L+M-1
          tmp = 0;
        if (n>=M)
          kmin = n-M+1;
        else
          kmin = 1;
        end

        if (n<L)
          kmax = n;
        else
          kmax = L;
        end
        %fprintf('kmin=%i,kmax=%i,n=%i\n',kmin,kmax,n);
        for k=kmin:kmax
          tmp = tmp + x(k)*imp_left(n-k+1);
        end
        x_conv(n)=tmp;
      end
  else
      x_conv = zeros(1,L+M-1);
     for i = 1:buffer_size:L
     
      for n=1:buffer_size+M-1
          tmp = 0;
        if ((n+i-1)>=M)
          kmin = n+i-1   -M+1;
        else
          kmin = 1;
        end

        if ((n+i-1)<L)
          kmax = n+i-1  ;
        else
          kmax = L;
        end
        %fprintf('kmin=%i,kmax=%i,n=%i\n',kmin,kmax,n);
        for k=kmin:kmax
          tmp = tmp + x(k)*imp_left(n-k+1);
        end
        x_conv(n+i-1)=tmp;
      end
     
     end
  end

  soundsc(x_conv,Fs);
else
  % FFT-base convolution (faster)
  if(overlap_add == 0)
      x_conv = fconv(x,imp_left);
      
  else
      x = [x' zeros(1,ceil(L/buffer_size)*buffer_size - L)];
      
      x_conv = [];
      temp=zeros(1,buffer_size+M-1);
      for i=1:buffer_size:L
          temp_overlap = [temp(buffer_size+1:end),zeros(1,buffer_size)];
          temp=fconv(x(i:(i+buffer_size-1)),imp_left') + temp_overlap;
          x_conv = [x_conv, temp(1:buffer_size)]; 
        
      end
      
  end
  soundsc(x_conv,Fs);
end



%% END