GNGlike = function(x0,n,pol,ce,V,Z, S ,y, s)
{
# Softmax likelihood for Go -No Go
 b = matrix(0,4,S);
 N = matrix(0,4,S);
 
  p[1]=x0[1];
  p[2]=x0[2];
  p[3]=x0[3];
  p[4]=x0[4];
  e     = exp(x0[5])/(1+exp(x0[5]));
  rho1  = exp(x0[6]);
  rho2  = exp(x0[7]);
  

 
 GnGlike=0;
 
    for( i in 1 :n)
    {
         t = V[i]; 
               
        N[2*t+1,Z[i]]= b[2*t+1,Z[i]];
        N[2*t+2,Z[i]]= b[2*t+2,Z[i]];
        
        T=sum(exp((N[(2*t+1):(2*t+2),Z[i]]+p[(2*t+1):(2*t+2)])));     
            # Calculate Probability of action k=1 
            # Collection
        GnGlike = GnGlike-log(exp(N[pol[i],Z[i]]+p[pol[i]]) /T);
    
        if( ce[i]>=0)
        {
            re = ce[i]*rho1;
        }
        else
        {
            re = ce[i]*rho2;
        }
        
        b[pol[i],Z[i]] = b[pol[i],Z[i]]+e*(re-b[pol[i],Z[i]]);  #Instrumental update
     
     
        
    }
   GnGlike = GnGlike -sum(-log(sqrt(2*pi*s[1:4]^2)) -(x0[1:4]-y[1:4])^2/(2*s[1:4]^2))-(-log(rho2*sqrt(2*pi*s[7]^2)) -(log(rho2)-y[7])^2/(2*s[7]^2));#Normal Priors
   GnGlike = GnGlike -(-log(e*sqrt(2*pi*s[5]^2)) -(log(e/(1-e))-y[5])^2/(2*s[5]^2)-log(1-e));
   GnGlike = GnGlike -(-log(rho1*sqrt(2*pi*s[6]^2)) -(log(rho1)-y[6])^2/(2*s[6]^2));
}

