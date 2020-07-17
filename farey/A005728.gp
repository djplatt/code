n1=1000000;
a=1;

printf("\nconst size_t A005728_table_size=%u;\n",1+n1);
printf("const unsigned long A005728_table[%u]={\n",1+n1);
for(n=1,n1,{
  printf("%uUL, // %u\n",a,n-1);
  a=a+eulerphi(n);
})

printf("%uUL, // %u\n};\n",a,n1);
quit;
