function [ Psi ] = FFT( psi,V,g,dt,k,x)
    m =1; hbar=1;
    dx = x(2)-x(1);dk = k(2)-k(1);
    M =length(x);Nx = 2*M;
    N = 1;% Define k-space grid

    
    psi=psi.*exp(-0.5*1i*dt*(V+(g/hbar)*abs(psi).^2));%coordinate space
    psi_k=fftshift(fft(psi)/Nx);
    psi_k=psi_k.*exp(-0.5*dt*1i*(hbar/m)*k.^2);%kinetic space
    psi=ifft(ifftshift(psi_k))*Nx;
    Psi=psi.*exp(-0.5*1i*dt*(V+(g/hbar)*abs(psi).^2));
    
end