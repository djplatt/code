do_k(lam,c0,c1,c2,c3,c)=(log(c3)-log((c1-2)*c2+c*log(2)/c0))/log(lam)+2;

\\ aiming for k=11

c0=0.6601618158;

\\ c1=7.8342; \\ H-B & Puchta
c1=7.8209;
c1=7.64;

c2=1.93656;
\\c2=1.93642;

\\ D=10000, mod 36465
\\c3=3.02846998;
\\ D=35537, mod 255255
\\c3=3.02857943365230;
\\wishful thinking
c3=3.02922591;

lam=0.8594000;


print("Unconditional with (9) k<=",do_k(lam,c0,c1,c2,c3,109/154));
print("Unconditional w/o (9)  k<=",do_k(lam,c0,c1,c2,c3,1));

\\ aiming for k=6
\\c3=3.01098859;
\\ D=14309 mod 255255
c3=3.01104895;
\\ wishful thinking
\\c3=3.01175656;
lam=0.7163436;

print("On GRH with (9) k<=",do_k(lam,c0,c1,c2,c3,1/2));
print("On GRH w/o (9)  k<=",do_k(lam,c0,c1,c2,c3,1));

