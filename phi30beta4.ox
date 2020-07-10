/**************************************************************************
 	 PROGRAM: phi30beta4.ox								                  
 																		  
 	 USAGE: Simulation of the test size -  
            corrected and uncorrected tests
			
 	 NULL HYPOTHESIS: H_0: \beta_4 = 0                                    								                              *
 																		  
 	 MODEL: g(mu,lambda) = X*betas   betas  = (beta_1,...,beta_p)	  
            fixed precision (\phi)                                    	  
 																		  
 	 AUTHOR: Cristine Rauber 									  																		  *
 **************************************************************************/

	// header files
	#include <oxstd.oxh>
	#include <oxprob.oxh>
	#import  <maximize>
	#import  <maxsqp>

	// global variables
	static decl y;		
	static decl N;
	static decl X;
	static decl Z;
	static decl Xr;
	static decl Xt;
	static decl Zt;
	static decl Xrt;
	static decl R=10000;   // Monte Carlo replications

	// irrestricted log-likelihood function 
	floglik(const vtheta, const adFunc, const avScore, const amHess)
	{

	decl beta    = vtheta[0:3];	  
	decl phi     = vtheta[4];
	decl lambda	 = vtheta[5];
	decl eta1	 = X*beta;
	decl mu		 = 1.0 - (1.0 + lambda .* exp(eta1)) .^ (-1.0 ./ lambda);
	decl p		 = mu .* phi;
	decl q		 = (1.0 - mu) .* phi;

	decl ystar   = log(y ./ (1.0 - y));
	decl ydag	 = log(1.0 - y);
	decl mustar	 = polygamma(mu .* phi, 0) - polygamma((1.0 - mu) .* phi, 0);
	decl mudag	 = polygamma((1.0 - mu) .* phi, 0) - polygamma(phi, 0);

	decl T       = diag( exp(eta1) .* (1.0 + lambda .* exp(eta1)) .^ (-(1.0 + (1.0 ./ lambda))) ); 
	decl H       = unit(N);	
	decl P       = phi .* unit(N);
	decl M       = diag(mu);
	decl rho     = (1.0 ./ lambda) .* ((1.0 ./ (exp(-eta1) + lambda)) -
                   (log(1.0 + lambda .* exp(eta1)) ./ lambda)) .* ((1.0 + lambda .* exp(eta1)) .^
				   (-1.0 ./ lambda));

	adFunc[0]    = double ( sumc( log(densbeta(y, p, q)) ) );

	// first order derivatives of the log-likelihood function
	if(avScore)
	{
	(avScore[0])[0:3] = Xt*P*T*(ystar-mustar);
	(avScore[0])[4]   = Zt*H*(M*(ystar-mustar)+(ydag-mudag));
	(avScore[0])[5]	  = rho'*P*(ystar-mustar);
	}

	return 1; // 1 indicates success 
	
	}

	// restricted log-likelihood function 	
	flogliknull(const vtheta, const adFunc, const avScore, const amHess)
	{

	decl beta     =	vtheta[0:2];   
	decl phi 	  =	vtheta[3];
	decl lambda	  = vtheta[4];
	decl eta1	  = Xr*beta;	  
	decl mu		  = 1.0 - (1.0 + lambda .* exp(eta1)) .^ (-1.0 ./ lambda);
	decl p		  = mu .* phi;
	decl q		  = (1.0 - mu) .* phi;

	decl ystar    = log(y ./ (1.0 - y));
	decl ydag	  = log(1.0 - y);
	decl mustar	  = polygamma(mu .* phi, 0) - polygamma((1.0 - mu) .* phi, 0);
	decl mudag	  = polygamma((1.0 - mu) .* phi, 0) - polygamma(phi, 0);

	decl T        = diag( exp(eta1) .* (1.0 + lambda .* exp(eta1)) .^ (-(1.0 + (1.0 ./ lambda))) ); 
	decl H        = unit(N);		  
	decl P        = phi .* unit(N);
	decl M        = diag(mu);
	decl rho      = (1.0 ./ lambda) .* ((1.0 ./ (exp(-eta1) + lambda)) -
	                (log(1.0 + lambda .* exp(eta1)) ./ lambda)) .* ((1.0 + lambda .* exp(eta1)) .^
					(-1.0 ./ lambda));

	adFunc[0]     = double ( sumc( log(densbeta(y, p, q)) ) );

	// first order derivatives of the log-likelihood function
	if(avScore)
	{
	(avScore[0])[0:2] = Xrt*P*T*(ystar-mustar);
	(avScore[0])[3]   = Zt*H*(M*(ystar-mustar)+(ydag-mudag));
	(avScore[0])[4]	  = rho'*P*(ystar-mustar);
	}
	
	return 1; // 1 indicates success 
	} 
	
	main()
	{

	// variables used in the maximization of the log-likelihood function
	decl dfunc0, dfunc1, conv0, conv1, vtheta0, vtheta1;
	decl vlo0, vlo1, vhi0, vhi1;

	// other variables used
	decl i, ir, j, fail, wneg, time, vEMVhat, vEMVtil, vstat;
	decl cv1, cv5, cv10, rej1w, rej5w, rej10w;
	decl rej1w1, rej5w1, rej10w1, rej1w2, rej5w2, rej10w2;
	decl theta, thetahat, thetatil, ystar, ydagger;

    // variables used in the model
	decl eta1, mu, phi, p, q, beta, lambda;

	oxwarning(0);
	time    = timer();     // time
	fail    = 0;		   // counter for the number of failures
	wneg    = 0;           // counter for the number of times when w<0
	vEMVtil = zeros(R,5);  // matrix to store the estimates under H0
	vEMVhat = zeros(R,6);  // matrix to store the estimates under H1
	vstat   = zeros(R,3);  // matrix to store the test statistics

	// significance levels
	cv1     = quanchi(0.99,1);
	cv5     = quanchi(0.95,1);
	cv10    = quanchi(0.9,1);

   	ranseed("MWC_52");	 // pseudorandom number generator
	ranseed(2018);       // seed of the pseudorandom sequence 

	for(j=20; j<=90; j+=10)	 // loop for sample sizes
	{

	N = j;	   // sample
	
	// generation of the model regressors
	X   = 1~(ranu(N,3));	 // uniform (0,1) 
	Z   = X[][0];	         // uniform (0,1) 
	Xr  = X[][0:2];		     // uniform (0,1) 
	Xt  = X';                // X transposed
	Zt  = Z';                // Z transposed
    Xrt = Xr';               // Xr transposed
	
	// parameters values
	beta    = <-1.5;1.5;4.0;0.0>;
	phi     = 30; 
	lambda	= 0.5;
	theta   = beta | phi | lambda;		// true parameter vector
	
    eta1  = X*beta;			 // linear predictor
	mu	  = 1.0 - (1.0 + lambda .* exp(eta1)) .^ (-1.0 ./ lambda);
	p	  =	mu .* phi;
	q	  =	(1.0 - mu) .* phi;
	
	for(ir=0; ir<R; ir++)		   // Monte Carlo loop
	{

	// sample generation
	y     = ranbeta(N, 1, p, q);

	// checking if y=1 or y=0
	for(i=0; i<N; i++)
	{
	if(y[i]==0.0000)
	   y[i]= 0.0001;
	else if(y[i]==1.0000)
	        y[i]= 0.9999;
	}
	
	ystar	 = log(y ./ (1.0 - y));
	ydagger	 = log(1.0 - y);

	// initial values
	vtheta1   = <-1.2;1.2;3.5;-3.5;9.5;0.2>;
	vtheta0   = <-1.2;1.2;3.5;9.5;0.2>;

	// boundaries for the initial values
	vlo1	= <-.Inf;-.Inf;-.Inf;-.Inf;0.001;0.001>;
	vhi1	= <+.Inf;+.Inf;+.Inf;+.Inf;+.Inf;10.000>;
	vlo0    = <-.Inf;-.Inf;-.Inf;0.001;0.001>;
	vhi0	= <+.Inf;+.Inf;+.Inf;+.Inf;10.000>;
	
	// convergence checking
	conv1 =	MaxSQP(floglik, &vtheta1, &dfunc1, 0, 0, 0, 0, vlo1, vhi1);
	conv0 =	MaxSQP(flogliknull, &vtheta0, &dfunc0, 0, 0, 0, 0, vlo0, vhi0);
	
	if( (conv1 == MAX_CONV  || conv1 == MAX_WEAK_CONV) && (conv0 == MAX_CONV || conv0 == MAX_WEAK_CONV) )
	{
	decl iota        = ones(N,1);		  // N-dimensional vector of ones
	decl Ystar       = diag(ystar);
	decl Ydagger	 = diag(ydagger);
	
	// quantities under H1***********************************************************************************************************
	decl eta1hat     = X*vtheta1[0:3];
	decl phihat      = vtheta1[4] ;
	decl lambdahat   = vtheta1[5];
	decl muhat       = 1.0 - (1.0 + lambdahat .* exp(eta1hat)) .^ (-1.0 ./ lambdahat);
	
	decl Hhat        = unit(N);
	decl That        = diag(exp(eta1hat) .* (1.0 + lambdahat .* exp(eta1hat)) .^ (-(1.0 + (1.0 ./ lambdahat))));
	decl Phat        = phihat .* unit(N);
	decl Muhat       = diag(muhat);
	decl Mustarhat	 = diag(polygamma(muhat .* phihat, 0) - polygamma((1.0 - muhat) .* phihat, 0)); // média de ystar
	decl Mudaggerhat = diag(polygamma((1.0 - muhat) .* phihat, 0) - polygamma(phihat, 0));	 // média de ydag
	decl Vstarhat    = diag(polygamma(muhat .* phihat, 1) + polygamma((1.0 - muhat) .* phihat, 1)); // covariância entre ystar e ystar
	decl Vdaggerhat  = diag(polygamma((1.0 - muhat) .* phihat, 1) - polygamma(phihat, 1));	// covariância entre ydag e ydag
	decl Chat        = diag(-polygamma((1.0 - muhat) .* phihat, 1)); // covariância entre ystar e ydag
	decl Shat        = diag((lambdahat - lambdahat .* (1.0 + lambdahat) .* (1.0 - muhat) .^ (lambdahat)) ./
	                       (((muhat - 1.0) .^ 2) .* ((1.0 - muhat) .^ (lambdahat) - 1.0) .^ 2)); // g''(mu, lambda)
	decl Qhat        = zeros(N);
	decl rhohat      = (1.0 ./ lambdahat) .* ((1.0 ./ (exp(-eta1hat) + lambdahat)) -
                       (log(1.0 + lambdahat .* exp(eta1hat)) ./ lambdahat)) .* ((1.0 + lambdahat .* exp(eta1hat)) .^
					   (-1.0 ./ lambdahat)); 
    decl varrhohat   = ((1.0 + lambdahat .* exp(eta1hat)) .^ (-2.0 - (1.0 ./ lambdahat)) .*
	                   (-exp(eta1hat) .* (lambdahat .^ 2) .* (2.0 + exp(eta1hat) .* (1.0 + 3.0 .* lambdahat)) +
					   (1.0 + lambdahat .* exp(eta1hat)) .* log(1.0 + lambdahat .* exp(eta1hat)) .*
					   (2.0 .* lambdahat .* (1.0 + exp(eta1hat) .* (1.0 + lambdahat)) -
					   (1.0 + lambdahat .* exp(eta1hat)) .* log(1.0 + lambdahat .* exp(eta1hat))))) ./
					   (lambdahat .^ 4);
	decl what        = (exp(eta1hat) .* (1.0 + lambdahat .* exp(eta1hat)) .^ (-2.0 - (1.0 ./ lambdahat)) .*
	                   (-exp(eta1hat) .* lambdahat .* (1.0 + lambdahat) + (1.0 + lambdahat .* exp(eta1hat)) .*
					   log(1.0 + lambdahat .* exp(eta1hat)))) ./ (lambdahat .^ 2);

	// observed information***********************************************************************************
	decl Jbbhat  =  Xt*(Phat*That*Vstarhat + That*Shat*That*(Ystar - Mustarhat))*That*Phat*X;
	decl Jbghat	 = -Xt*((Ystar - Mustarhat) - Phat*(Muhat*Vstarhat + Chat))*That*Hhat*Z;
	decl Jblhat	 =  Xt*(Phat*Vstarhat*Phat*That*rhohat - Phat*(Ystar - Mustarhat)*what);
	decl Jgbhat	 =  Jbghat';
	decl Jgghat	 =  Zt*(Hhat*(Muhat*Vstarhat*Muhat + (Muhat + Muhat)*Chat + Vdaggerhat)+
	                   (Muhat*(Ystar - Mustarhat) + (Ydagger - Mudaggerhat))*Hhat*Qhat*Hhat)*Hhat*Z;
	decl Jglhat	 = -Zt*((Ystar - Mustarhat) - Phat*(Muhat*Vstarhat + Chat))*Hhat*rhohat;
	decl Jlbhat	 =  Jblhat';
	decl Jlghat	 =  Jglhat';
	decl Jllhat  =  ((Phat .^ 2)*Vstarhat*(rhohat .^ 2) - Phat*(Ystar - Mustarhat)*varrhohat)'*iota;

	decl Jhat    =  (Jbbhat~Jbghat~Jblhat) | (Jgbhat~Jgghat~Jglhat) | (Jlbhat~Jlghat~Jllhat);
	decl invJhat =  invert(Jhat);  // inverse of Jhat

	// Fisher's information***********************************************************************************
	decl Kbbhat  =  Xt*Phat*That*Vstarhat*That*Phat*X;
	decl Kbghat	 =	Xt*Phat*(Muhat*Vstarhat + Chat)*Hhat*That*Z;
	decl Kblhat	 =	Xt*Phat*Vstarhat*Phat*That*rhohat;
	decl Kgbhat	 =	Kbghat';
	decl Kgghat	 =	Zt*Hhat*(Muhat*Vstarhat*Muhat + (Muhat + Muhat)*Chat + Vdaggerhat)*Hhat*Z;
	decl Kglhat	 =	Zt*Phat*(Muhat*Vstarhat + Chat)*Hhat*rhohat;
	decl Klbhat	 =	Kblhat';
	decl Klghat	 =	Kglhat';
	decl Kllhat	 =	rhohat'*(Phat .^ 2)*Vstarhat*rhohat;
			 
	decl Khat    =  (Kbbhat~Kbghat~Kblhat) | (Kgbhat~Kgghat~Kglhat) | (Klbhat~Klghat~Kllhat);
	decl invKhat =  invert(Khat); // inversa de Khat
	
	// quantities under H0***********************************************************************************
	decl eta1til     = X*(vtheta0[0:2] | theta[3]);
	decl phitil      = vtheta0[3];
	decl lambdatil   = vtheta0[4];
	decl mutil       = 1.0 - (1.0 + lambdatil .* exp(eta1til)) .^ (-1.0 ./ lambdatil);

	decl Htil        = unit(N);
	decl Ttil        = diag(exp(eta1til) .* (1.0 + lambdatil .* exp(eta1til)) .^ (-(1.0 + (1.0 ./ lambdatil))));
	decl Ptil        = phitil .* unit(N);
	decl Mutil       = diag(mutil);
	decl mustartil   = polygamma(mutil .* phitil, 0) - polygamma((1.0 - mutil) .* phitil, 0);
	decl Mustartil	 = diag(mustartil); // média de ystar
	decl mudaggertil = polygamma((1.0 - mutil) .* phitil, 0) - polygamma(phitil, 0);
	decl Mudaggertil = diag(mudaggertil);	 // média de ydag
	decl Vstartil    = diag(polygamma(mutil .* phitil, 1) + polygamma((1.0 - mutil) .* phitil, 1)); // covariância entre ystar e ystar
	decl Vdaggertil  = diag(polygamma((1.0 - mutil) .* phitil, 1) - polygamma(phitil, 1));	// covariância entre ydag e ydag
	decl Ctil        = diag(-polygamma((1.0 - mutil) .* phitil, 1)); // covariância entre ystar e ydag
	decl Stil        = diag((lambdatil - lambdatil .* (1.0 + lambdatil) .* (1.0 - mutil) .^ (lambdatil)) ./
	                       (((mutil - 1.0) .^ 2) .* ((1.0 - mutil) .^ (lambdatil) - 1.0) .^ 2)); // g''(mu, lambda)
	decl Qtil        = zeros(N);
	decl rhotil      = (1.0 ./ lambdatil) .* ((1.0 ./ (exp(-eta1til) + lambdatil)) -
                       (log(1.0 + lambdatil .* exp(eta1til)) ./ lambdatil)) .* ((1.0 + lambdatil .* exp(eta1til)) .^
					   (-1.0 ./ lambdatil)); 
    decl varrhotil   = ((1.0 + lambdatil .* exp(eta1til)) .^ (-2.0 - (1.0 ./ lambdatil)) .*
	                   (-exp(eta1til) .* (lambdatil .^ 2) .* (2.0 + exp(eta1til) .* (1.0 + 3.0 .* lambdatil)) +
					   (1.0 + lambdatil .* exp(eta1til)) .* log(1.0 + lambdatil .* exp(eta1til)) .*
					   (2.0 .* lambdatil .* (1.0 + exp(eta1til) .* (1.0 + lambdatil)) -
					   (1.0 + lambdatil .* exp(eta1til)) .* log(1.0 + lambdatil .* exp(eta1til))))) ./
					   (lambdatil .^ 4);
	decl wtil        = (exp(eta1til) .* (1.0 + lambdatil .* exp(eta1til)) .^ (-2.0 - (1.0 ./ lambdatil)) .*
	                   (-exp(eta1til) .* lambdatil .* (1.0 + lambdatil) + (1.0 + lambdatil .* exp(eta1til)) .*
					   log(1.0 + lambdatil .* exp(eta1til)))) ./ (lambdatil .^ 2);
					   
    // observed information***********************************************************************************
	decl Jbbtil  =  Xt*(Ptil*Ttil*Vstartil + Ttil*Stil*Ttil*(Ystar - Mustartil))*Ttil*Ptil*X;
	decl Jbgtil	 = -Xt*((Ystar - Mustartil) - Ptil*(Mutil*Vstartil + Ctil))*Ttil*Htil*Z;
	decl Jbltil	 =  Xt*(Ptil*Vstartil*Ptil*Ttil*rhotil - Ptil*(Ystar - Mustartil)*wtil);
	decl Jgbtil	 =  Jbgtil';
	decl Jggtil	 =  Zt*(Htil*(Mutil*Vstartil*Mutil + (Mutil + Mutil)*Ctil + Vdaggertil)+
	                   (Mutil*(Ystar - Mustartil) + (Ydagger - Mudaggertil))*Htil*Qtil*Htil)*Htil*Z;
	decl Jgltil	 = -Zt*((Ystar - Mustartil) - Ptil*(Mutil*Vstartil + Ctil))*Htil*rhotil;
	decl Jlbtil	 =  Jbltil';
	decl Jlgtil	 =  Jgltil';
	decl Jlltil  =  ((Ptil .^ 2)*Vstartil*(rhotil .^ 2) - Ptil*(Ystar - Mustartil)*varrhotil)'*iota;

	decl Jtil    =  (Jbbtil~Jbgtil~Jbltil) | (Jgbtil~Jggtil~Jgltil) | (Jlbtil~Jlgtil~Jlltil);
	decl invJtil =  invert(Jtil); // inverse of Jtil

	// Fisher's information***********************************************************************************
	decl Kbbtil  =  Xt*Ptil*Ttil*Vstartil*Ttil*Ptil*X;
	decl Kbgtil	 =	Xt*Ptil*(Mutil*Vstartil + Ctil)*Ttil*Htil*Z;
	decl Kbltil	 =	Xt*Ptil*Vstartil*Ptil*Ttil*rhotil;
	decl Kgbtil	 =	Kbgtil';
	decl Kggtil	 =	Zt*Htil*(Mutil*Vstartil*Mutil + (Mutil + Mutil)*Ctil + Vdaggertil)*Htil*Z;
	decl Kgltil	 =	Zt*Ptil*(Mutil*Vstartil + Ctil)*Htil*rhotil;
	decl Klbtil	 =	Kbltil';
	decl Klgtil	 =	Kgltil';
	decl Klltil	 =	rhotil'*(Ptil .^ 2)*Vstartil*rhotil;
			 
	decl Ktil    =  (Kbbtil~Kbgtil~Kbltil) | (Kgbtil~Kggtil~Kgltil) | (Klbtil~Klgtil~Klltil);
	decl invKtil =  invert(Ktil);  // inverse of Ktil

	// score function under H0***********************************************************************************
	decl escorebetatil	  =	 Xt*Ptil*Ttil*(ystar - mustartil);
	decl escoregamatil	  =	 Zt*Htil*(Mutil*(ystar - mustartil) + (ydagger - mudaggertil));
	decl escorelambdatil  =	 rhotil'*Ptil*(ystar - mustartil);

	decl escoretil        =  escorebetatil | escoregamatil | escorelambdatil;
	
	// qbar***********************************************************************************
	decl qbeta	 =	Xt*Phat*That*(Vstarhat*(Phat*Muhat - Ptil*Mutil) + (Phat - Ptil)*Chat)*iota;
	decl qgama	 =	Zt*Hhat*((Muhat*Vstarhat + Chat)*(Phat*Muhat - Ptil*Mutil) +
	                        (Muhat*Chat + Vdaggerhat)*(Phat - Ptil))*iota;
	decl qlambda =	rhohat'*Phat*(Vstarhat*(Phat*Muhat - Ptil*Mutil) + Chat*(Phat - Ptil))*iota;

	decl qbar    =  qbeta | qgama | qlambda;
	
	// upsilonbar***********************************************************************************
	decl upsbb   =	Xt*Phat*That*Vstarhat*Ttil*Ptil*X;
	decl upsbg	 =	Xt*Phat*That*(Vstarhat*Mutil + Chat)*Htil*Z;
	decl upsbl	 =	Xt*Phat*That*Vstarhat*Ptil*rhotil;
	decl upsgb	 =	Zt*Hhat*(Muhat*Vstarhat + Chat)*Ttil*Ptil*X;
	decl upsgg	 =	Zt*Hhat*(Muhat*Vstarhat*Mutil + (Muhat + Mutil)*Chat + Vdaggerhat)*Htil*Z;
	decl upsgl	 =	Zt*Hhat*(Muhat*Vstarhat + Chat)*Ptil*rhotil;
	decl upslb	 =	rhohat'*Phat*Vstarhat*Ttil*Ptil*X;
	decl upslg	 =	rhohat'*Phat*(Vstarhat*Mutil + Chat)*Htil*Z;
	decl upsll	 =	rhohat'*Phat*Vstarhat*Ptil*rhotil;
	
	decl upsbar	    = (upsbb~upsbg~upsbl) | (upsgb~upsgg~upsgl) | (upslb~upslg~upsll);
	decl invupsbar	= invert(upsbar); // inverse of upsbar
	
	decl nuiJtil    = (Jtil[0:2][0:2] ~ Jtil[0:2][4:5]) | (Jtil[4:5][0:2] ~ Jtil[4:5][4:5]);
	decl nuisance   = Ktil*invupsbar*Jhat*invKhat*upsbar;
	decl nuisance2  = (nuisance[0:2][0:2] ~ nuisance[0:2][4:5]) | (nuisance[4:5][0:2] ~ nuisance[4:5][4:5]);

	// Likelihood ratio test statistics***********************************************************************************

	vstat[ir][0]     = 2*(dfunc1-dfunc0); // likelihood ratio test statistic
	
	// checking if w<0
	if(vstat[ir][0] <= 0)
	{
	wneg++;
	ir--;
	}
   
	decl epson      = fabs( ((fabs(determinant(Ktil))*fabs(determinant(Khat))*fabs(determinant(nuiJtil)))^(0.5))/
						   (fabs(determinant(upsbar))*fabs(determinant(nuisance2))^(0.5))*
						   (fabs(escoretil'*invupsbar*Khat*invJhat*upsbar*invKtil*escoretil)^(1/2))/
						   fabs((vstat[ir][0])^((1/2)-1.0)*escoretil'*invupsbar*qbar) );

	// checking if w<=0.1					   
	if(vstat[ir][0] <= 0.1)
	{

	vstat[ir][1]     = vstat[ir][0];	                    
	vstat[ir][2]     = vstat[ir][0];
	
	}
    else
	{

	vstat[ir][1]     = vstat[ir][0]-2*log(epson);  // Skovgaard's modified likelihood ratio test statistic w*
	vstat[ir][2]     = vstat[ir][0]*(1.0-log(epson)/vstat[ir][0])^2; // Skovgaard's modified likelihood ratio test statistic w*

	}

	vEMVtil[ir][]    = vtheta0';	 // parameter estimates under H0
	vEMVhat[ir][]    = vtheta1';	 // parameter estimates under H1

	}
	else
	{
	fail++;
	ir--;
	}

	} // end of Monte Carlo loop

	// null rejection rates
    rej1w   = (sumc(vstat[][0] .> cv1)/R)*100;
	rej5w   = (sumc(vstat[][0] .> cv5)/R)*100;
	rej10w  = (sumc(vstat[][0] .> cv10)/R)*100;
	
	rej1w1  = (sumc(vstat[][1] .> cv1)/R)*100;
	rej5w1  = (sumc(vstat[][1] .> cv5)/R)*100;
	rej10w1 = (sumc(vstat[][1] .> cv10)/R)*100;
	
	rej1w2  = (sumc(vstat[][2] .> cv1)/R)*100;
	rej5w2  = (sumc(vstat[][2] .> cv5)/R)*100;
	rej10w2 = (sumc(vstat[][2] .> cv10)/R)*100;

	thetahat = (meanc(vEMVhat))';
	thetatil = (meanc(vEMVtil))';
	thetatil = thetatil[0:2] | zeros(1,1) | thetatil[3:4];

	// save the test statistics to a .mat file
	if(N==20)
	{
	savemat("phi30n20b4.mat", vstat);
	}
	else if(N==30)
	{
	savemat("phi30n30b4.mat", vstat);
	}
	else if(N==40)
	{
	savemat("phi30n40b4.mat", vstat);
	}
	else if(N==50)
	{
	savemat("phi30n50b4.mat", vstat);
	}
	else if(N==60)
	{
	savemat("phi30n60b4.mat", vstat);
	}
	else if(N==70)
	{
	savemat("phi30n70b4.mat", vstat);
	}
	else if(N==80)
	{
	savemat("phi30n80b4.mat", vstat);
	}
	else if(N==90)
	{
	savemat("phi30n90b4.mat", vstat);
	}
	//**************************************************************************************************************************
	
	// printing results
	println("------------------------------RESULTS----------------------------");
	println("\t\t\t Program:", oxfilename(0));
	println("\t\t\t OX version:", oxversion());
	println("\t\t\t Pseudorandom number generator:", ranseed(""));
	println("\t\t\t Seed: 2018");
	println("\t\t\t Sample size:", N);
	println("\t\t\t Monte Carlo replications:", R);
	println("\t\t\t Monte Carlo replications with failures:", fail);
	println("\t\t\t Monte Carlo replications with w<=0:", wneg);
	println("\t\t\t Optimization algorithm applied: MaxSQP");
	println("%10.4f", "%c", {"theta","hat", "til"}, theta~thetahat~thetatil);
	println("-----------------------------------------------------------------");
	println("\t\t\t\t NULL REJECTION RATES");
	println("-----------------------------------------------------------------");
	println("%12.2f", "%c", {"1%", "5%", "10%"}, "%r", {"w", "w*", "w**"},
			 (rej1w~rej5w~rej10w) | (rej1w1~rej5w1~rej10w1) | (rej1w2~rej5w2~rej10w2));
	println("-----------------------------------------------------------------");

	}  // end of loop for sample sizes
	
	println("\t\t\t Date:", date());
	println("\t\t\t Time:", timespan(time));
	println("-----------------------------------------------------------------");
				
	}  // end of main****************************************************************