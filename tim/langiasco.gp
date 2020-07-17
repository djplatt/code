/*** 	A. LANGUASCO & A. ZACCAGNINI  ***/
/***  	Implementation of Pintz-Ruzsa's method paper  ***/
/***  	(as described in their paper paper on Acta Arith. 109,(2003))   ***/

/***  Improved by K. BELABAS using the gexpo function ***/


\\ global variables

global(L, L1, majpsiU, minpsiU, majpsiL, minpsiL, momentmajU, momentminU, momentmajL, momentminL);



\\ Goal function and its first derivative : gam e gamprime
\\ gamma is internally used for the Euler Gamma function

{gam(x,lambda) =
exp(lambda*cos(x));
}


{gamprime(x,lambda) = 
(-1)*lambda*sin(x)*exp(lambda*cos(x));
}


\\ gexpo is a libpari function;  returns the binary exponent of x or 
\\ the maximal binary exponent of the coefficients of x.
\\ lg is a libpari function;  returns the length of x

install(gexpo, lG);


\\ main function: it realizes the approximation via the matrix method
\\ c is the level 
\\ k is the degree of the polynomials
\\ lambda is the point in which we evaluate psi
\\ numdigits is the  precision for the final result

{PintzRuzsa(c, k, lambda, numdigits) = 

\\ local variables

local(log2 = log(2));

local(c1 = c*log2);

local(eU, eL, expU, expL, eps2);

local(sigmaupper, sigmalower, sigma0upper,sigma0lower, rhoupper, delta);

local(rholower, Uupper, Ulower, Bupper, Blower, K, j, l, eps1);

local(i, i1, i2, m, n0, j1, j2, alphaupp, alphalow, betalow,betaupp); 


\\ Building the coefficients of the polynomials

K = k+1;

sigmaupper = vector(k);
sigmalower = vector(k);
rhoupper = vector(k);
rholower = vector(k);

\\ print(sigmaupper);
\\ print(sigmalower);
\\ print(rhoupper);
\\ print(rholower);

sigma0upper = 0;
sigma0lower = 0;

for(i=1,K, 
		i1=2*Pi*(i-1)/K; i2=i1+Pi/K;
		sigma0upper += gam(i1,lambda);
    	sigma0lower += gam(i2,lambda);
    	);
    	    	  	
sigma0upper /= K;
sigma0lower /= K;

for (j=1, k,
 	for(i=1,K, 
 		i1=2*Pi*(i-1)/K; i2=i1+Pi/K;
 		sigmaupper[j] +=  gam(i1,lambda) * cos(j*i1);
 		sigmalower[j] +=  gam(i2,lambda) * cos(j*i2);
 		rhoupper[j] +=  gamprime(i1,lambda) * sin(j*i1);
 		rholower[j] +=  gamprime(i2,lambda) * sin(j*i2)
 		);
 	sigmaupper[j] /= K;
 	sigmalower[j] /= K;
 	rhoupper[j] /= K;
 	rholower[j] /= K
 	);
 		
\\print(lambda);
\\print(sigma0upper);
\\print(sigma0lower);
\\print(sigmaupper-sigmalower);
\\print(sigmalower);
\\error("");
\\print(rhoupper);
\\print(rholower);
\\error("");

Bupper = vector(k+1);
Blower = vector(k+1);

Bupper[1] = sigma0upper;
Blower[1] = sigma0lower;


for (i=1, k,
 		Bupper[i+1]= (2/K) * ((K-i)*sigmaupper[i]-rhoupper[i]);
 		Blower[i+1]= (2/K) * ((K-i)*sigmalower[i]-rholower[i]);
 	);

\\ print("Bupper =", Bupper);
\\ printp("Blower =", Blower);


\\ Building the matricies U

m = k-1;
Uupper = matrix(m+1,m+1);
Ulower = matrix(m+1,m+1);

\\ print("k= ", k); 

\\ This definition follows from line -6 and -5 at p.178 pf Pintz-Ruzsa paper.

for (j=1,m+1,
	n0=(j-1)%2;
	forstep (n=n0,m, 2,
			j1=(j-1+n)\2; j2= abs(j-1-n)\2; 
			Uupper[n+1,j1+1] += Bupper[j]/2;
			Uupper[n+1,j2+1] += Bupper[j]/2;
			Ulower[n+1,j1+1] += Blower[j]/2;
			Ulower[n+1,j2+1] += Blower[j]/2;
		)
	);	
	
\\print("lambda=",lambda);
\\print("Uupper =", Uupper);
\\print("Ulower =", Ulower);


\\ eps1 is the desired precision

eps1= 10^(-(numdigits));

\\ print(eps1);

L=0;
delta=10;

\\ expU and expL are used to keep track of the exponents used

expU = 0;
expL = 0;
eps2 = eps1^3;

until (delta<eps1,

		\\ first normalization step:  reducing the size of the exponents;
		\\ if the maximal exponent eU (eL) in the matrices is too large
		\\ we reduce it by reducing every element of 10^eU (10^eL)
		\\ via shifting of eU (el) position and we keep track of it
		\\ increasing the global exponent expU (exp L)
		
        eU = gexpo(Uupper); if (eU, Uupper >>= eU; expU += eU);
        eL = gexpo(Ulower); if (eL, Ulower >>= eL; expL += eL);

		\\ second normalization step: the entry in the lower right corner
		\\ decreases too quickly comparing to the others
		
        if (abs(Uupper[m+1,m+1]) < eps2, Uupper[m+1,m+1] = 0);
        if (abs(Ulower[m+1,m+1]) < eps2, Ulower[m+1,m+1] = 0);
	
	\\ squaring out the matrices and keeping track of the exponents

	Uupper = Uupper^2; expU *= 2;
	Ulower = Ulower^2; expL *= 2;
	
	\\ keeping track of the number of squares
	
	L++;
    
    \\ approximations, see Lemma 7, p. 179, of Pintz-Ruzsa paper

\\if(L==10,print("L=",L);print("U^L=",Uupper));
    
	alphaupp = Uupper[1,1];
	betaupp = sum(i=1,m+1, Uupper[1,i]);
	
\\if(L==10,print("alphaupp = ", alphaupp);print("betaupp = ", betaupp));
	\\ computing the precision
	
	delta = log(betaupp / alphaupp) / (lambda*2^L);
	
	\\ print("delta= ", delta);
	);

\\ computing the  moments and the approximations
L1 = 2^L;

\\print("expU=",expU);
\\print("L1=",L1);
\\print("L=",L);
\\print("lambda=",lambda);
\\print("betaupp=",betaupp);

momentmajU = (log(betaupp) + expU * log2)/L1;
momentminU = (log(alphaupp)+ expU * log2)/L1;

\\print("momentmajU=",momentmajU);
alphalow = Ulower[1,1];
betalow = sum(i=1,m+1, Ulower[1,i]);
	
momentmajL = (log(betalow) + expL * log2)/L1;
momentminL = (log(alphalow)+ expL * log2)/L1;

\\ computing the approximations for the function psi of Pintz-Ruzsa paper

majpsiU= (momentmajU + c1)/lambda;
minpsiU= (momentminU + c1)/lambda;
majpsiL= (momentmajL + c1)/lambda;
minpsiL= (momentminL + c1)/lambda;


\\ if the degree of the polynomials is too small we are not able to reach
\\ the fixed precision

if (abs(majpsiU-minpsiL) > eps1, 
			error("increase the degree of the polynomials")
			);

\\ return the approximated value of psi(lambda)

return(majpsiU)
}


\\ This function realizes a dyadic approximation using the previous function
\\ it  starts using the approximated values of psi in 1-T, 1, 1+T
\\ and it moves to a new set of three points to look for a saddle point
\\ which is approximated using a precision 10^(-numdigits)

\\ c is the level 
\\ k is the degree of the polynomials
\\ numdigits is the  precision for the final result


{PintzRuzsa_psiapprox(c,k,numdigits) = local(iteration = 0,eps0,r0,r1,r2,l0,T);
	
	numdigits += 2;
	
	default(realprecision,numdigits+10);        \\ set working precision;
	
	iteration=0;
	
	print("The expected number of correct decimal digits is = ", numdigits-2);
	
	eps0=10^(-numdigits); 
	
	print("Approximation for the minimal lambda is = 10^(-", numdigits-2,")");
	
	r0=0; r1=0; r2=0;
	
	T=1/10;
	l0=1;
	
	until (T<eps0,	 
		iteration++;
		\\ print("l0-T = ", l0-T);
		r0 = PintzRuzsa(c,k,l0-T,numdigits);
	    \\ print ("r0 = ", r0);
		\\ print("l0 = ", l0);		 
		r1 = PintzRuzsa(c,k,l0,numdigits); 
		\\ print ("r1 = ", r1);
		\\ print("l0+T = ", l0+T);
		r2 = PintzRuzsa(c,k,l0+T,numdigits);
	    \\ print ("r2 = ", r2);
		until ((r0>r1), iteration++;
					l0=l0-2*T; 
					r2=r1; 
					r1=r0; 
					\\ print("l0 = ", l0); 
					r0 = PintzRuzsa(c,k,l0,numdigits);			
				    \\ print ("r0 = ", r0);
				);
		until ((r2>r1),iteration++;
					l0=l0+2*T; 
					r0=r1; 
					r1=r2; 
					\\ print("l0 = ", l0); 
					r2 = PintzRuzsa(c,k,l0,numdigits);
				    \\ print ("r2 = ", r2);
				);
		T /= 10;
		);
		
	\\ output of the final results	
		
	print ("Needed matrix exponent for this precision is = 2^", L);
	print ("Number of iterations in the dyadic procedure = ", iteration);
	print ("The approximated final values are: ");
	print ("max-upper-matrix = ", majpsiU);
	print ("min-upper-matrix = ", minpsiU);
	print ("max-lower-matrix = ", majpsiL);
	print ("min-lower-matrix = ", minpsiL);
	print ("The approximated values for the moments are: ");
	print ("moment-max-upper-matrix = ", momentmajU);
	print ("moment-min-upper-matrix = ", momentminU);
	print ("moment-max-lower-matrix = ", momentmajL);
	print ("moment-min-lower-matrix = ", momentminL);
	print ("The minimal lambda is in [", (l0-T)*1.0, ",",(l0+T)*1.0,"]");
	print ("Final result (in the centre of the last interval): d = ", r1);
	}


/*******************************************************

HOW TO RUN THIS PROGRAM

- with gp2c (faster) to be run as 



1) > gp2c-run -pmy_ -g -W /path/PRmethodfinal.gp 
2) (10:01) gp > PintzRuzsa_psiapprox(1/2,14,12)

- with gp (slowler) to be run as 
1) gp
2) \r /path/PRmethodfinal.gp
3) (10:01) gp > PintzRuzsa_psiapprox(1/2,14,12)

WHERE  PATH IS THE path of the file 
you saved PRmethodfinal.gp

******************************************************

Results:
-----------------------------------------------------
4/5 on a quad core pc
-----------------------------------------------------
[languasc@labsrv0 ~]$ nice /usr/local/Gruppi/PariGP/bin/gp2c-run -pmy_ -g -W PRmethodfinal.gp 
PRmethodfinal.gp.c: In function `my_PintzRuzsa':
PRmethodfinal.gp.c:94: warning: unused variable `my_j'
PRmethodfinal.gp.c:94: warning: unused variable `my_l'
PRmethodfinal.gp.c:94: warning: unused variable `my_i'
                     GP/PARI CALCULATOR Version 2.3.2 (released)
            amd64 running linux (x86-64/GMP-4.2.2 kernel) 64-bit version
           compiled: Nov 30 2007, gcc-3.4.3 20041212 (Red Hat 3.4.3-9.EL4)
                  (readline v4.3 enabled, extended help available)

                       Copyright (C) 2000-2006 The PARI Group

PARI/GP is free software, covered by the GNU General Public License, and comes 
WITHOUT ANY WARRANTY WHATSOEVER.

Type ? for help, \q to quit.
Type ?12 for how to get moral (and possibly technical) support.

parisize = 8000000, primelimit = 500000

-----------------------------------------------------
10 digits
-----------------------------------------------------

? PintzRuzsa_psiapprox(4/5+10^(-20),13,10)
The expected number of correct decimal digits is = 10
Approximation for the minimal lambda is = 10^(-10)
Needed matrix exponent for this precision is = 2^39
Number of iterations in the dyadic procedure = 162
The approximated final values are: 
max-upper-matrix = 0.9123781030527322323638
min-upper-matrix = 0.9123781030521349899098
max-lower-matrix = 0.9123781030527020894117
min-lower-matrix = 0.9123781030521048469576
The approximated values for the moments are: 
moment-max-upper-matrix = 0.6181690752736390706245
moment-min-upper-matrix = 0.6181690752728714301565
moment-max-lower-matrix = 0.6181690752736003276488
moment-min-lower-matrix = 0.6181690752728326871808
The minimal lambda is in [1.285307939545900000000,1.285307939546100000000]
Final result (in the centre of the last interval): d = 0.9123781030527322323638
time = 5,609 ms.
 

-----------------------------------------------------
20 digits
-----------------------------------------------------

? PintzRuzsa_psiapprox(4/5+10^(-20),20,20)
The expected number of correct decimal digits is = 20
Approximation for the minimal lambda is = 10^(-20)
Needed matrix exponent for this precision is = 2^72
Number of iterations in the dyadic procedure = 280
The approximated final values are: 
max-upper-matrix = 0.91237810305275834972134056121053
min-upper-matrix = 0.91237810305275834972127103303804
max-lower-matrix = 0.91237810305275834972133828300078
min-lower-matrix = 0.91237810305275834972126875482828
The approximated values for the moments are: 
moment-max-upper-matrix = 0.61816907526618862709860869006595
moment-min-upper-matrix = 0.61816907526618862709851932495382
moment-max-lower-matrix = 0.61816907526618862709860576186487
moment-min-lower-matrix = 0.61816907526618862709851639675273
The minimal lambda is in [1.2853079395377972468345900000000,1.2853079395377972468346100000000]
Final result (in the centre of the last interval): d = 0.91237810305275834972134056121053
time = 49,239 ms.
 

-----------------------------------------------------
30 digits
-----------------------------------------------------

?  PintzRuzsa_psiapprox(4/5+10^(-20),27,30)
The expected number of correct decimal digits is = 30
Approximation for the minimal lambda is = 10^(-30)
Needed matrix exponent for this precision is = 2^105
Number of iterations in the dyadic procedure = 424
The approximated final values are: 
max-upper-matrix = 0.912378103052758349721358590459910416545890
min-upper-matrix = 0.912378103052758349721358590459902322401313
max-lower-matrix = 0.912378103052758349721358590459910399454390
min-lower-matrix = 0.912378103052758349721358590459902305309813
The approximated values for the moments are: 
moment-max-upper-matrix = 0.618169075266188626931299358124411230167273
moment-min-upper-matrix = 0.618169075266188626931299358124400826698984
moment-max-lower-matrix = 0.618169075266188626931299358124411208199432
moment-min-lower-matrix = 0.618169075266188626931299358124400804731143
The minimal lambda is in [1.28530793953779724665119741228003900000000,1.28530793953779724665119741228004100000000]
Final result (in the centre of the last interval): d = 0.912378103052758349721358590459910416545890
time = 4mn, 19,948 ms.
 

-----------------------------------------------------
50 digits
-----------------------------------------------------

? PintzRuzsa_psiapprox(4/5+10^(-20),39,50)
The expected number of correct decimal digits is = 50
Approximation for the minimal lambda is = 10^(-50)
Needed matrix exponent for this precision is = 2^172
Number of iterations in the dyadic procedure = 624
The approximated final values are: 
max-upper-matrix = 0.91237810305275834972135859045990929414085615244675930311308271
min-upper-matrix = 0.91237810305275834972135859045990929414085615244675924826502561
max-lower-matrix = 0.91237810305275834972135859045990929414085615244675930308170857
min-lower-matrix = 0.91237810305275834972135859045990929414085615244675924823365148
The approximated values for the moments are: 
moment-max-upper-matrix = 0.61816907526618862693129935808253892939804025384544671244136441
moment-min-upper-matrix = 0.61816907526618862693129935808253892939804025384544664194472116
moment-max-lower-matrix = 0.61816907526618862693129935808253892939804025384544671240103899
moment-min-lower-matrix = 0.61816907526618862693129935808253892939804025384544664190439573
The minimal lambda is in [1.2853079395377972466511974122341479975582119999426043900000000,1.2853079395377972466511974122341479975582119999426044100000000]
Final result (in the centre of the last interval): d = 0.91237810305275834972135859045990929414085615244675930311308271
? ##
  ***   last result computed in 30mn, 26,314 ms.
 
-----------------------------------------------------

-----------------------------------------------------
2/3 on a quad core pc
-----------------------------------------------------
[languasc@labsrv0 ~]$ nice /usr/local/Gruppi/PariGP/bin/gp2c-run -pmy_ -g -W PRmethodfinal.gp 
PRmethodfinal.gp.c: In function `my_PintzRuzsa':
PRmethodfinal.gp.c:94: warning: unused variable `my_j'
PRmethodfinal.gp.c:94: warning: unused variable `my_l'
PRmethodfinal.gp.c:94: warning: unused variable `my_i'
                     GP/PARI CALCULATOR Version 2.3.2 (released)
            amd64 running linux (x86-64/GMP-4.2.2 kernel) 64-bit version
           compiled: Nov 30 2007, gcc-3.4.3 20041212 (Red Hat 3.4.3-9.EL4)
                  (readline v4.3 enabled, extended help available)

                       Copyright (C) 2000-2006 The PARI Group

PARI/GP is free software, covered by the GNU General Public License, and comes 
WITHOUT ANY WARRANTY WHATSOEVER.

Type ? for help, \q to quit.
Type ?12 for how to get moral (and possibly technical) support.

parisize = 8000000, primelimit = 500000

-----------------------------------------------------
10 digits
-----------------------------------------------------

? PintzRuzsa_psiapprox(2/3+10^(-20),12,10)
The expected number of correct decimal digits is = 10
Approximation for the minimal lambda is = 10^(-10)
Needed matrix exponent for this precision is = 2^39
Number of iterations in the dyadic procedure = 160
The approximated final values are: 
max-upper-matrix = 0.8337213168426473838898
min-upper-matrix = 0.8337213168420918764748
max-lower-matrix = 0.8337213168425383672048
min-lower-matrix = 0.8337213168419828597898
The approximated values for the moments are: 
moment-max-upper-matrix = 0.4413015839592436256564
moment-min-upper-matrix = 0.4413015839586416916237
moment-max-lower-matrix = 0.4413015839591254978852
moment-min-lower-matrix = 0.4413015839585235638525
The minimal lambda is in [1.083575154049900000000,1.083575154050100000000]
Final result (in the centre of the last interval): d = 0.8337213168426473838898
time = 4,673 ms.
 

-----------------------------------------------------
20 digits
-----------------------------------------------------

? PintzRuzsa_psiapprox(2/3+10^(-20),19,20)
The expected number of correct decimal digits is = 20
Approximation for the minimal lambda is = 10^(-20)
Needed matrix exponent for this precision is = 2^72
Number of iterations in the dyadic procedure = 277
The approximated final values are: 
max-upper-matrix = 0.83372131684338485515459227435255
min-upper-matrix = 0.83372131684338485515452760477806
max-lower-matrix = 0.83372131684338485515458851235650
min-lower-matrix = 0.83372131684338485515452384278202
The approximated values for the moments are: 
moment-max-upper-matrix = 0.44130158394251203319483400833145
moment-min-upper-matrix = 0.44130158394251203319476393398731
moment-max-lower-matrix = 0.44130158394251203319482993192600
moment-min-lower-matrix = 0.44130158394251203319475985758187
The minimal lambda is in [1.0835751540289729521761900000000,1.0835751540289729521762100000000]
Final result (in the centre of the last interval): d = 0.83372131684338485515459227435255
time = 41,941 ms.
 

-----------------------------------------------------
30 digits
-----------------------------------------------------

? PintzRuzsa_psiapprox(2/3+10^(-20),26,30)
The expected number of correct decimal digits is = 30
Approximation for the minimal lambda is = 10^(-30)
Needed matrix exponent for this precision is = 2^105
Number of iterations in the dyadic procedure = 423
The approximated final values are: 
max-upper-matrix = 0.833721316843384855154592152588205415890590
min-upper-matrix = 0.833721316843384855154592152588197887361263
max-lower-matrix = 0.833721316843384855154592152588205404464508
min-lower-matrix = 0.833721316843384855154592152588197875935182
The approximated values for the moments are: 
moment-max-upper-matrix = 0.441301583942512033013203743652502539297678
moment-min-upper-matrix = 0.441301583942512033013203743652494381570354
moment-max-lower-matrix = 0.441301583942512033013203743652502526916660
moment-min-lower-matrix = 0.441301583942512033013203743652494369189336
The minimal lambda is in [1.08357515402897295195834526956583900000000,1.08357515402897295195834526956584100000000]
Final result (in the centre of the last interval): d = 0.833721316843384855154592152588205415890590
time = 3mn, 53,681 ms.


-----------------------------------------------------
50 digits
-----------------------------------------------------

? PintzPintzRuzsa_psiapprox(2/3+10^(-20),37,50)
The expected number of correct decimal digits is = 50
Approximation for the minimal lambda is = 10^(-50)
Needed matrix exponent for this precision is = 2^172
Number of iterations in the dyadic procedure = 628
The approximated final values are: 
max-upper-matrix = 0.83372131684338485515459215258820372013640387733149037455160858
min-upper-matrix = 0.83372131684338485515459215258820372013640387733149032353630951
max-lower-matrix = 0.83372131684338485515459215258820372013640387733149037427262730
min-lower-matrix = 0.83372131684338485515459215258820372013640387733149032325732823
The approximated values for the moments are: 
moment-max-upper-matrix = 0.44130158394251203301320374362026095967995267147401959779762169
moment-min-upper-matrix = 0.44130158394251203301320374362026095967995267147401954251871114
moment-max-lower-matrix = 0.44130158394251203301320374362026095967995267147401959749532451
moment-min-lower-matrix = 0.44130158394251203301320374362026095967995267147401954221641397
The minimal lambda is in [1.0835751540289729519583452695271703132722006062074439900000000,1.0835751540289729519583452695271703132722006062074440100000000]
Final result (in the centre of the last interval): d = 0.83372131684338485515459215258820372013640387733149037455160858
time = 26mn, 25,158 ms.
? 



********************************************/

