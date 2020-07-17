g++ -O2 -static -march=nocona -msse3 -finline-functions -fomit-frame-pointer check_proth.cpp -I ~/cln-1.3.2/include/ -L ~/cln-1.3.2/src/.libs -lcln -o check_proth
