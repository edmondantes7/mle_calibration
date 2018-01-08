function [ theta, sigma, mu, loglhd ] = CalibrateOrnsteinUhlenbeckMaxLikelihood(S, deltat) 
% Calibrate an OU process by maximum likelihood. 
  
%% Reference 
% Based on the algorithm and software described at : 
% http://www.sitmo.com/doc/Calibrating_the_Ornstein-Uhlenbeck_model 
  n = length(S)-1; 
  
  Sx  = sum( S(1:end-1) ); 
  Sy  = sum( S(2:end) ); 
  Sxx = sum( S(1:end-1).^2 ); 
  Sxy = sum( S(1:end-1).*S(2:end) ); 
  Syy = sum( S(2:end).^2 ); 
  
  theta  = (Sy*Sxx - Sx*Sxy) / ( n*(Sxx - Sxy) - (Sx^2 - Sx*Sy) ); 
   
  mu = -(1/deltat)*log((Sxy - theta*Sx - theta*Sy + n*theta^2) / (Sxx -2*theta*Sx + n*theta^2)); 
  alpha  = exp(-  mu*deltat); 
  alpha2 = exp(-2*mu*deltat); 
  sigmahat2 = (1/n)*(Syy - 2*alpha*Sxy + alpha2*Sxx - ... 
               2*theta*(1-alpha)*(Sy - alpha*Sx) + n*theta^2*(1-alpha)^2); 
  sigma = sqrt(sigmahat2*2*mu/(1-alpha2)); 
  aux1 = S(2:end) - S(1:end-1)*exp(-mu*deltat)- theta*(1-exp(-mu*deltat));
  aux2 = aux1.^2;
  loglhd = (-n/2*log(2*pi)-n/2*log(sigmahat2) -(1/(2*sigmahat2))*sum(aux2))/n;
  theta
  sigma
  mu
  loglhd
end 