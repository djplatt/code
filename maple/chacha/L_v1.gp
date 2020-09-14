eps=0;
/*

Integrate L''/L'(s) around the box 0->1/2+eps->1/2+eps+t0*I->t0*I->0
where t0=max(3,6/log(q))
Divide the imaginary part of the answer by 2Pi and round to an integer

This should give N-P where N is number of zeros of L' in the box and P is
the number of poles of L', counted with multiplicity/order resp.

We know that L' for primitive chi is entire, so P=0

*/

/* if we call lfun(.,small,.) we get an error */
mylfun(L,s,d)=if(abs(s)<1e-10,lfun(L,0,d),lfun(L,s,d));

chacha(q,n)={local(L,t0,r,rr);t0=max(6/log(q),3);L=lfuninit(Mod(n,q),[1/4+eps/2,1/4+eps/2,2*t0],2);rr=(intnum(t=0,1/2+eps,imag(mylfun(L,t-t0*I,2)/mylfun(L,t-t0*I,1)-mylfun(L,t+t0*I,2)/mylfun(L,t+t0*I,1)))+intnum(t=-t0,t0,real(mylfun(L,1/2+eps+I*t,2)/mylfun(L,1/2+eps+I*t,1)-mylfun(L,I*t,2)/mylfun(L,I*t,1))))/(2*Pi);r=round(rr);if(abs(r-rr)>0.01,printf("I=%f ",rr));printf("q=%d: chi=%d: N-P=%d\n",q,n,r)}

/* do all primitive characters q in [3,22] 
for(q=3,22,if(q%4!=2,G=znstar(q,1);for(chi=2,q-1,if(gcd(q,chi)==1,if(zncharconductor(G,chi)==q,chacha(q,chi))))));

/* do all primitive odd characters q in [23,215]  
for(q=23,215,if(q%4!=2,G=znstar(q,1);for(chi=2,q-1,if(gcd(q,chi)==1,if(zncharconductor(G,chi)==q,if(zncharisodd(G,chi),chacha(q,chi)))))));
*/

/* here are the problems (i.e. not within 0.01 of an int) */

allocatemem(1000000000);
\p 100
probs=[163,162,167,70,15,14,149,93,181,31,197,56,191,55,179,69,127,29,173,18,181,95,163,50,209,17,197,5,151,15,109,33,181,72,199,74,13,5,101,63,203,45,199,149,193,140,157,22,179,86,13,2,211,133,199,134,19,7];

/*for(n=1,length(probs)/2,chacha(probs[2*n-1],probs[2*n]));*/

/* try higher precision */

\p 500
probs1=[149,93,19,7,163,162,15,14,199,134,157,22,167,70];
for(n=1,length(probs1)/2,chacha(probs1[2*n-1],probs1[2*n]));


quit;

/* all go away at 500 bits, except Mod(10,11) which has a zero at s=0.17380922... */


/* show that Mod(10,11) has only the one zero inside gamma, the real one */
eps=0.1;
chacha1(q,n)={local(L,t0,r,rr);t0=max(6/log(q),3);L=lfuninit(Mod(n,q),[1/4,1/4,2*t0],2);rr=(intnum(t=0,1/2,imag(mylfun(L,t-eps*I,2)/mylfun(L,t-eps*I,1)-mylfun(L,t+t0*I,2)/mylfun(L,t+t0*I,1)))+intnum(t=-eps,t0,real(mylfun(L,1/2+I*t,2)/mylfun(L,1/2+I*t,1)-mylfun(L,I*t,2)/mylfun(L,I*t,1))))/(2*Pi);r=round(rr);if(abs(r-rr)>0.01,printf("I=%f ",rr));printf("q=%d: chi=%d: N-P=%d\n",q,n,r)}

chacha1(11,10);


/* find primitive order 2 (real) with L'(s)=0 in [0,1/2] */
for(q=3,300,if(q%4!=2,G=znstar(q,1);for(n=2,q-1,if(gcd(n,q)==1,if(znorder(Mod(n,q))==2,if(zncharconductor(G,n)==q,L=lfuninit(Mod(n,q),[1/4,1/4,6],1);if(lfun(L,0,1)*lfun(L,1/2,1)<0,printf("%d %d %f\n",q,n,solve(t=0,1/2,mylfun(L,t,1))))))))));



quit;
