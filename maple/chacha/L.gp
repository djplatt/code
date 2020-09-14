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

chacha(q,n)={local(L,t0,r,rr);t0=max(6/log(q),3);L=lfuninit(Mod(n,q),[1/4+eps/2,1/4+eps/2,2*t0],2);rr=(intnum(t=0,1/2+eps,imag(mylfun(L,t,2)/mylfun(L,t,1)-mylfun(L,t+t0*I,2)/mylfun(L,t+t0*I,1)))+intnum(t=0,t0,real(mylfun(L,1/2+eps+I*t,2)/mylfun(L,1/2+eps+I*t,1)-mylfun(L,I*t,2)/mylfun(L,I*t,1))))/(2*Pi);r=round(rr);if(abs(r-rr)>0.01,printf("I=%f ",rr));printf("q=%d: chi=%d: N-P=%d\n",q,n,r)}

/* do all primitive characters q in [3,215] */
for(q=3,215,if(q%4!=2,G=znstar(q,1);for(chi=2,q-1,if(gcd(q,chi)==1,if(zncharconductor(G,chi)==q,chacha(q,chi))))));
quit;
