
// assumes |x|<=1e-5 and returns log(1+x)
int_double log1p_small(int_double x)
{
  static bool init=false;
  static int_double err,third;
  if(!init)
    {
      init=true;
      third=1.0;
      third/=3.0;
      err=1.0;
      err/=5e25; // 1/5|x|^5
      err.left=err.right;
    }
  int_double x2,x3,x4,res;
  x2=sqr(x);
  x4=sqr(x2);
  x3=x2*x;
  res=err-x4*0.25;
  res+=x3*third;
  res-=x2*0.5;
  return res+x;
}
