from mpmath import *
import sys
mp.dps = 40
fn=sys.argv[1]
st=int(sys.argv[2])
en=int(sys.argv[3])
f=open(fn,'w')
for n in range(st,en):
   x=zetazero(n)
   w=zeta(x,1,1)
   y=im(x)
   z=str(y)
   print(n,z)
   f.write(z)
   f.write(' ')
   f.write(str(re(w)))
   f.write(' ')
   f.write(str(im(w)))
   f.write('\n')
f.close()
