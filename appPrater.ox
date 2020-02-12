/**************************************************************************
 	 PROGRAM: appPrater.ox								                  
 																		  
 	 USAGE: Computation of the modified test statistics   
            for the gasoline yield data from Prater (1956)
			
 	 NULL HYPOTHESIS: H_0: lambda = 1                                    								                              *
 																		  
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
	static decl Xrr;
	static decl Xrrt;
	static decl Z;
	static decl Xt;
	static decl Zt;

	// irrestricted log-likelihood function 
	floglik(const vtheta, const adFunc, const avScore, const amHess)
	{

	decl r       = columns(X);
	decl beta    = vtheta[0:(r-1)];	  
	decl phi     = vtheta[r];
	decl lambda	 = vtheta[r+1];
	decl eta1	 = X*beta;		  
	decl mu		 = 1.0 - (1.0 + lambda .* exp(eta1)) .^ (-1.0 ./ lambda);
	decl p		 = mu .* phi;
	decl q		 = (1.0 - mu) .* phi;

	decl ystar   = log(y ./ (1.0 - y));
	decl ydag	 = log(1.0 - y);
	decl mustar	 = polygamma(mu .* phi, 0) - polygamma((1.0 - mu) .* phi, 0);
	decl mudag	 = polygamma((1.0 - mu) .* phi, 0) - polygamma(phi, 0);

	decl T       = diag( exp(eta1) .* (1.0 + lambda .* exp(eta1)) .^ (-(1.0 + (1.0 ./ lambda))) ); // 1/g'(mu, lambda)
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
	(avScore[0])[0:(r-1)]  = Xt*P*T*(ystar-mustar);
	(avScore[0])[r]        = Zt*H*(M*(ystar-mustar)+(ydag-mudag));
	(avScore[0])[r+1]      = rho'*P*(ystar-mustar);
	}

    if( isnan(adFunc[0]) || isdotinf(adFunc[0]) ) 
 	return 0; 
    else
    return 1; // 1 indicates success 	
	
	}

	// restricted log-likelihood function 
	flogliknull(const vtheta, const adFunc, const avScore, const amHess)
    {
	
	decl r       = columns(X);
	decl beta    = vtheta[0:(r-1)];	  
	decl phi     = vtheta[r];
	decl lambda  = 1;
	decl eta1	 = X*beta;            
	decl mu	     = 1.0 - (1.0 + lambda .* exp(eta1)) .^ (-1.0 ./ lambda);
	decl p		 = mu .* phi;
	decl q		 = (1.0 - mu) .* phi;

	decl ystar   = log(y ./ (1.0 - y));
	decl ydag	 = log(1.0 - y);
	decl mustar  = polygamma(mu .* phi, 0) - polygamma((1.0 - mu) .* phi, 0);
	decl mudag	 = polygamma((1.0 - mu) .* phi, 0) - polygamma(phi, 0);

	decl T       = diag( exp(eta1) .* (1.0 + lambda .* exp(eta1)) .^ (-(1.0 + (1.0 ./ lambda))) ); // 1/g'(mu, lambda)
 	decl H       = unit(N);	
	decl P       = phi .* unit(N);
	decl M       = diag(mu);
	
    adFunc[0]    = double ( sumc( log(densbeta(y, p, q)) ) );

	// first order derivatives of the log-likelihood function
	if(avScore)
	{
	(avScore[0])[0:(r-1)]  = Xt*P*T*(ystar-mustar);
	(avScore[0])[r]        = Zt*H*(M*(ystar-mustar)+(ydag-mudag));
	}

	if( isnan(adFunc[0]) || isdotinf(adFunc[0]) ) 
 	return 0; 
    else
    return 1; // 1 indicates success 	
    }

	// log-likelihood function of the null model
	flogliknullaranda(const vtheta, const adFunc, const avScore, const amHess)
    {
	
	 decl beta    = vtheta[0];	  
	 decl lambda  = 1;		  
	 decl phi     = vtheta[1];
	 decl eta1	  = Xrr*beta;      
	 decl mu	  = 1.0 - (1.0 + lambda .* exp(eta1)) .^ (-1.0 ./ lambda);
	 decl p		  = mu .* phi;
	 decl q		  = (1.0 - mu) .* phi;

	 decl ystar   = log(y ./ (1.0 - y));
	 decl ydag	  = log(1.0 - y);
	 decl mustar  = polygamma(mu .* phi, 0) - polygamma((1.0 - mu) .* phi, 0);
	 decl mudag	  = polygamma((1.0 - mu) .* phi, 0) - polygamma(phi, 0);
	
	 decl T       = diag( exp(eta1) .* (1.0 + lambda .* exp(eta1)) .^ (-(1.0 + (1.0 ./ lambda))) ); // 1/g'(mu, lambda)

     adFunc[0]    = double ( sumc( log(densbeta(y, p, q)) ) );

	 // first order derivatives of the log-likelihood function
     if(avScore)
     {  
          (avScore[0])[0] = phi*Xrrt*T*(ystar-mustar); 
	      (avScore[0])[1] = double( sumc( mu .* (ystar-mustar) + (ydag-mudag) ) ); 
     }	     	  
 
     if( isnan(adFunc[0]) || isdotinf(adFunc[0]) ) 
 	 return 0; 
     else
     return 1; // 1 indicates success 	
    }

	main()
	{

	// variables used in the maximization of the log-likelihood function
	decl dfunc0, dfunc1, dfuncr, conv0, conv1, conv2, dscore;
	decl vtheta0, vtheta1, vthetar;
	decl vlo1, vhi1, vlo0, vhi0, vlo2, vhi2;

	// other variables used
	decl ybar, yvar, ystar, ydagger;
	decl r, s, k, gl, pseudoR2LR;
	decl w, pvw, ws, pvws, wss, pvwss;

	// variables used for the initial values
	decl betaols, gamaols, phiols, varols, muols, etaols, lambdaini;

	// variables used in the model
	decl data, batch, temp;

	oxwarning(0);
	data     = loadmat("gasoline.mat");	 // load the data
	y        = data[][10];			     // variable of interest
	temp     = data[][9];				 // covariate temp
	batch    = data[][0:8];				 // covariate batch
	X        = 1~batch~temp;	         // matrix 32x11
	Z        = X[][0];				     // matrix 32x1 of 1's
	Xt       = X';                       // X transposed
	Zt       = Z';                       // Z transposed

	k  = 1;                              // number of parameters of interest
	r  = columns(X);				     // number of parameters in the mean submodel
	s  = columns(Z);				     // number of parameters in the precision submodel
	N  = rows(data);			         // sample size
	gl = N-(r+s);					     // degrees of fredom

	ystar	 = log(y ./ (1.0 - y));	     // transformed variable
	ydagger	 = log(1.0 - y);			 // transformed variable

	ols2c(ystar, X, &betaols);           // store the ols estimates in betaols
	etaols    = X*betaols;
	muols     = exp(etaols) ./ (1.0 + exp(etaols));
	varols    = ((ystar - etaols)' * (ystar - etaols)) ./ ((N - r) * ((1 ./ (muols .* (1.0 - muols))) .^ (2)));
	phiols   = double( meanc((muols .* (1.0 - muols) ./ varols) - 1.0) );
	lambdaini = 1;	  // initial value for lambda (logit)

	// initial values
	vtheta1   = betaols | phiols  | lambdaini;
	vtheta0   = betaols | phiols;

	// boundaries for the initial values
	vlo1    = <-.Inf;-.Inf;-.Inf;-.Inf;-.Inf;-.Inf;-.Inf;-.Inf;-.Inf;-.Inf;-.Inf;0.001;0.001>;
	vhi1	= <+.Inf;+.Inf;+.Inf;+.Inf;+.Inf;+.Inf;+.Inf;+.Inf;+.Inf;+.Inf;+.Inf;+.Inf;+.Inf>;
	vlo0    = <-.Inf;-.Inf;-.Inf;-.Inf;-.Inf;-.Inf;-.Inf;-.Inf;-.Inf;-.Inf;-.Inf;0.001>;
	vhi0	= <+.Inf;+.Inf;+.Inf;+.Inf;+.Inf;+.Inf;+.Inf;+.Inf;+.Inf;+.Inf;+.Inf;+.Inf>;

	ybar = meanc(y);   // mean of y
    yvar = varc(y);	   // variance of y
	
	println("-----------------------------------------------------------------");
	println("\t\t\t\t BETA REGRESSION ESTIMATION");
	println("-----------------------------------------------------------------");
     
	println("\n MEAN AND VARIANCE OF Y:\n ", "%10.5f", ybar~yvar);  
    println("\n INITIAL VALUES FOR THE ML ESTIMATION:\n ", "%16.5f", vtheta1);
	println("-----------------------------------------------------------------");

	// convergence checking
	conv1 =	MaxSQP(floglik, &vtheta1, &dfunc1, 0, 0, 0, 0, vlo1, vhi1);
	conv0 =	MaxSQP(flogliknull, &vtheta0, &dfunc0, 0, 0, 0, 0, vlo0, vhi0);
	println("\n CONVERGENCE STATUS UNDER H1: ", MaxConvergenceMsg(conv1));
	println("\n CONVERGENCE STATUS UNDER H0: ", MaxConvergenceMsg(conv0));
 
	if( (conv1 == MAX_CONV  || conv1 == MAX_WEAK_CONV) && (conv0 == MAX_CONV || conv0 == MAX_WEAK_CONV) )
	{
	
	decl iota        = ones(N,1);  // N-dimensional vector of ones
	decl Ystar       = diag(ystar);
	decl Ydagger	 = diag(ydagger);
	decl r           = columns(X);
	
	// quantities under H1***********************************************************************************************************
	decl eta1hat     = X*vtheta1[0:(r-1)];
	decl phihat      = vtheta1[r];
	decl lambdahat   = vtheta1[r+1];
	decl muhat       = 1.0 - (1.0 + lambdahat .* exp(eta1hat)) .^ (-1.0 ./ lambdahat);
	
	decl Hhat        = unit(N);
	decl That        = diag(exp(eta1hat) .* (1.0 + lambdahat .* exp(eta1hat)) .^ (-(1.0 + (1.0 ./ lambdahat))));
	decl Phat        = phihat .* unit(N);
	decl Muhat       = diag(muhat);
	decl mustarhat   = polygamma(muhat .* phihat, 0) - polygamma((1.0 - muhat) .* phihat, 0);
	decl Mustarhat	 = diag(polygamma(muhat .* phihat, 0) - polygamma((1.0 - muhat) .* phihat, 0)); 
	decl Mudaggerhat = diag(polygamma((1.0 - muhat) .* phihat, 0) - polygamma(phihat, 0));
	decl vstarhat    = polygamma(muhat .* phihat, 1) + polygamma((1.0 - muhat) .* phihat, 1);
	decl Vstarhat    = diag(polygamma(muhat .* phihat, 1) + polygamma((1.0 - muhat) .* phihat, 1)); 
	decl Vdaggerhat  = diag(polygamma((1.0 - muhat) .* phihat, 1) - polygamma(phihat, 1));	
	decl Chat        = diag(-polygamma((1.0 - muhat) .* phihat, 1)); 
	decl Shat        = diag((lambdahat - lambdahat .* (1.0 + lambdahat) .* (1.0 - muhat) .^ (lambdahat)) ./
	                      (((muhat - 1.0) .^ 2) .* ((1.0 - muhat) .^ (lambdahat) - 1.0) .^ 2)); 
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
	decl Jbbhat  =  Xt*(Phat*That*Vstarhat + Shat*(That^2)*(Ystar - Mustarhat))*That*Phat*X;
	decl Jbghat	 = -Xt*((Ystar - Mustarhat) - Phat*(Muhat*Vstarhat + Chat))*That*Hhat*Z;
	decl Jblhat	 =  Xt*((Phat^2)*Vstarhat*That*rhohat - Phat*(Ystar - Mustarhat)*what);
	decl Jgbhat	 =  Jbghat';
	decl Jgghat	 =  Zt*(Hhat*(Muhat*Vstarhat*Muhat + (Muhat + Muhat)*Chat + Vdaggerhat)+
	                   (Muhat*(Ystar - Mustarhat) + (Ydagger - Mudaggerhat))*(Hhat^2)*Qhat)*Hhat*Z;
	decl Jglhat	 = -Zt*((Ystar - Mustarhat) - Phat*(Muhat*Vstarhat + Chat))*Hhat*rhohat;
	decl Jlbhat	 =  Jblhat';
	decl Jlghat	 =  Jglhat';
	decl Jllhat  =  ((Phat^2)*Vstarhat*(rhohat .^ 2) - Phat*(Ystar - Mustarhat)*varrhohat)'*iota;

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
	decl Kllhat	 =	rhohat'*(Phat^2)*Vstarhat*rhohat;
			 
	decl Khat    =  (Kbbhat~Kbghat~Kblhat) | (Kgbhat~Kgghat~Kglhat) | (Klbhat~Klghat~Kllhat);
	decl invKhat =  invert(Khat); // inverse of Khat

	// quantities under H0***********************************************************************************
	decl eta1til     = X*vtheta0[0:(r-1)];
	decl phitil      = vtheta0[r];
	decl lambdatil   = 1;
	decl mutil       = 1.0 - (1.0 + lambdatil .* exp(eta1til)) .^ (-1.0 ./ lambdatil);
	
	decl Htil        = unit(N);
	decl Ttil        = diag(exp(eta1til) .* (1.0 + lambdatil .* exp(eta1til)) .^ (-(1.0 + (1.0 ./ lambdatil))));
	decl Ptil        = phitil .* unit(N);
	decl Mutil       = diag(mutil);
	decl mustartil   = polygamma(mutil .* phitil, 0) - polygamma((1.0 - mutil) .* phitil, 0);
	decl Mustartil	 = diag(mustartil); 
	decl mudaggertil = polygamma((1.0 - mutil) .* phitil, 0) - polygamma(phitil, 0);
	decl Mudaggertil = diag(mudaggertil);	
	decl Vstartil    = diag(polygamma(mutil .* phitil, 1) + polygamma((1.0 - mutil) .* phitil, 1)); 
	decl Vdaggertil  = diag(polygamma((1.0 - mutil) .* phitil, 1) - polygamma(phitil, 1));
	decl Ctil        = diag(-polygamma((1.0 - mutil) .* phitil, 1)); 
	decl Stil        = diag((lambdatil - lambdatil .* (1.0 + lambdatil) .* (1.0 - mutil) .^ (lambdatil)) ./
	                       (((mutil - 1.0) .^ 2) .* ((1.0 - mutil) .^ (lambdatil) - 1.0) .^ 2)); 
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
	decl invJtil =  invert(Jtil); // inverse  of Jtil

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
	
	decl nuiJtil    = Jtil[0:(r+s-1)][0:(r+s-1)];
	decl nuisance   = Ktil*invupsbar*Jhat*invKhat*upsbar;
	decl nuisance2  = nuisance[0:(r+s-1)][0:(r+s-1)];

	// Likelihood ratio test statistics***********************************************************************************

	          w     = 2*(dfunc1-dfunc0); 	 // likelihood ratio test statistic

             pvw    = 1.0-probchi(w,1);		 // p-value of w

	decl epson      = fabs( ((fabs(determinant(Ktil))*fabs(determinant(Khat))*fabs(determinant(nuiJtil)))^(0.5))/
						   (fabs(determinant(upsbar))*fabs(determinant(nuisance2))^(0.5))*
						   (fabs(escoretil'*invupsbar*Khat*invJhat*upsbar*invKtil*escoretil)^(1/2))/
						   fabs((w)^((1/2)-1.0)*escoretil'*invupsbar*qbar) );
	
	 ws			  = w-2*log(epson);		     // Skovgaard's modified likelihood ratio test statistic w*
	 pvws		  = 1.0-probchi(ws,1);		 // p-value of w*
	 wss		  = w*(1.0-log(epson)/w)^2;	 // Skovgaard's modified likelihood ratio test statistic w**
	 pvwss		  = 1.0-probchi(wss,1);		 // p-value of w**

	 // pseudoR2 based on likelihood******************************************************************************************
    Xrr           = ones(N, 1);
	Xrrt          = Xrr';
	decl ystarbar = meanc(ystar);
	decl lambdar  = 1;
	decl muhatr   = 1.0 - (1.0 + lambdar .* exp(ystarbar)) .^ (-1.0 ./ lambdar);

	// initial values (constant mean and precision, fixed lambda)***************************************************
	vthetar       = ystarbar | ((1.0/(varc(ystar)*muhatr*(1.0-muhatr))));
	vlo2          = <-.Inf;-.Inf>;
	vhi2		  =	<+.Inf;+.Inf>;

	// convergence checking
    conv2         = MaxSQP(flogliknullaranda, &vthetar, &dfuncr, 0, 0, 0, 0, vlo2, vhi2);                                                           
	                                                                                                           
    if(conv2 == MAX_CONV || conv2 == MAX_WEAK_CONV)                                                                    
    {                                                                                                       
      pseudoR2LR  = 1.0 - (exp(dfuncr)/exp(dfunc1))^(2/N);                                                
    }
	   	
	// measures of quality of the fitted model****************************************************************************************************
	decl pseudoR2  = (correlation(eta1hat~ystar)[0][1])^2;
	decl AIC       = -2*dfunc1+2*(r+s+1);
	decl BIC       = -2*dfunc1+(r+s+1)*log(N);
	//**************************************************************************************************************************

	// printing results
	println("\n PARAMETER ESTIMATES AND ASYMPTOTIC STANDARD ERRORS: ");
	decl stderrors = sqrt(diagonal(invKhat))';  // standard erros
	decl zstats    = vtheta1 ./ stderrors; 		// z test statistic
	println("%16.5f", "%c", {"estimates", "std. errors", "z stats", "p-values"}, "%r",
	                        {"intercept", "batch1", "batch2", "batch3", "batch4", "batch5", "batch6", "batch7",
							 "batch8", "batch9", "temp", "phi", "lambda"}, 
	                          vtheta1~stderrors~zstats~2.0*(1.0-probn(fabs(zstats))));
	println("\t Sample size:", N);	// sample size
	println("%r", {"pseudoR2", "pseudoR2LR", "AIC", "BIC"},	   // fitted model quality measures
	                pseudoR2  | pseudoR2LR | AIC | BIC);
	println("\n ASYMPTOTIC COVARIANCE MATRIX OF ML ESTIMATES:"); 
	println("%14.5f", invKhat);

	println("-----------------------------------------------------------------");
	println("\t\t\t\t LIKELIHOOD RATIO TEST STATISTICS");
	println("-----------------------------------------------------------------");
	println("\t\t\t NULL HYPOTHESIS: lambda = 1 (logit link):");
	println("-----------------------------------------------------------------");
	println("%16.6f", "%c", {"\t test statistic", "p-value"}, "%r", {"w", "w*", "w**"},
	                                            (w | ws | wss)~(pvw | pvws | pvwss));

	}
	else
	{
	println("\n\n ERROR: NO CONVERGENCE!\n\n");
	}

	println("-----------------------------------------------------------------");
	println("\t\t\t Program:", oxfilename(0));
	println("\t\t\t OX version:", oxversion());
	println("\t\t\t Optimization algorithm applied: MaxSQP");
	println("\t\t\t Date:", date());
	println("\t\t\t Time: ", time());
	println("-----------------------------------------------------------------");
				
	}  // end of main****************************************************************************************************

