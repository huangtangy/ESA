%% Solves the 1D GPE in imaginary time propagation
function [psi,mu] = get_ground_state(psi,dt,g1d,x,m,V,N0)

hbar = 1;%1.054e-34;
dx = x(2)-x(1); 
Nx=length(x);
dk=2*pi/(Nx*dx); M =Nx/2;
k =(-Nx/2:Nx/2-1)*dk;

j=1;
mu_error =1;
mu_old = 1;
psi_mid_old = psi(Nx/2); 

while abs(mu_error) > 1e-6
    
    psi = psi.*exp(-0.5*dt*(V+(g1d/hbar)*abs(psi).^2)); %Apply PE operator
    psi_k = fftshift(fft(psi))/Nx;
    psi_k = psi_k.*exp(-0.5*dt*(hbar/m)*k.^2); %Apply KE operator
    psi = ifft(ifftshift(psi_k))*Nx;
    psi = psi.*exp(-0.5*dt*(V+(g1d/hbar)*abs(psi).^2)); %Apply PE operator again
    
    
    psi_mid = psi((Nx)/2);
    mu = -log(psi_mid_old/psi_mid)/dt; %g1d*n0((Nx)/2)+D_0;%
    mu_error = abs((mu-mu_old)/mu); %Defines relative change in chemical potential
    
    %mu_error_old = mu_error;
    
    psi_mid_old = psi(Nx/2); 
    mu_old = mu; 
    
    if mod(j,50000) == 0
        mu_error
    end
    if j > 1e8 %Catches divergences
        'no solution found'
        break
    end
    psi = psi*sqrt(N0)/sqrt((dx*norm(psi).^2)); %Renormalises to keep atom number const
    
    %{
    plot(x,abs(psi).^2,'b-')
    drawnow
    %}

    j=j+1;
end 

end
