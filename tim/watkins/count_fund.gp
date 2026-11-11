count_fundamental(N) = {
  my(tot_1 = 0, tot_0 = 0, lim = sqrtint(N));
  my(mu = vector(lim, i, moebius(i)));
  
  \\ Case 1: D = 1 mod 4 -> n = 3 mod 4 (up to N)
  for(d = 1, lim,
    if(mu[d] == 0 || d % 2 == 0, next);
    my(max_m = N \ (d * d));
    if(max_m >= 3,
      tot_1 += mu[d] * ((max_m - 3) \ 4 + 1);
    );
  );
  
  \\ Case 2: D = 0 mod 4 -> n = 4k where k = 1 or 2 mod 4 (k up to N\4)
  my(N4 = N \ 4);
  my(lim4 = sqrtint(N4));
  for(d = 1, lim4,
    if(moebius(d) == 0 || d % 2 == 0, next);
    my(max_m = N4 \ (d * d));
    \\ Count m = 1 mod 4
    if(max_m >= 1, tot_0 += moebius(d) * ((max_m - 1) \ 4 + 1));
    \\ Count m = 2 mod 4
    if(max_m >= 2, tot_0 += moebius(d) * ((max_m - 2) \ 4 + 1));
  );
  
  print("Discriminants == 1 mod 4: ", tot_1);
  print("Discriminants == 0 mod 4: ", tot_0);
  print("Total Fundamental Discriminants: ", tot_1 + tot_0);
}
