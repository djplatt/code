first_zero(soln)={local(pos);pos=1;while(soln[pos]>0,pos++);pos};
try(order)={local(ts,tn,fz);ts=soln;for(n=1,5,tn=order[n];fz=first_zero(ts);if(fz>14-tn-1,return(0));if(ts[fz+tn+1]>0,return(0));ts[fz]=tn;ts[fz+tn+1]=tn);return(ts);}

/* put the 1 in pos 1 */
soln=[1,0,1,0,0,0,0,0,2,0,0,2,0,0];
forperm(5,p,q=p;for(pos=1,5,q[pos]+=2);res=try(q);if(res,print(res)));

/* try it in pos 5 */
soln=[0,0,1,0,1,0,0,0,2,0,0,2,0,0];
forperm(5,p,q=p;for(pos=1,5,q[pos]+=2);res=try(q);if(res,print(res)));


