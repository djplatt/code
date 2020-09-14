x=10^22;
w=55*10^14;
A=x-w;
B=x;

powers(n,A,B)=
{
   local (low,high);

   low=ceil(log(A)/n);
   high=floor(log(B)/n);
   for(i=ceil(exp(log(A)/n)),floor(exp(log(B)/n)),if(isprime(i),print(i," ",n)));
};

main()=
{

   local (max_n);

   max_n=floor(log(B)/log(2));

  for(i=2,max_n,powers(i,A,B));
}
   