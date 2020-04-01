# Feather
Feather: Lightweight Multi-party Updatable Delegated Private Set Intersection

# Dependencies

* GMP: https://gmplib.org/
* Cryptopp: https://www.cryptopp.com
* NTL: https://www.shoup.net/ntl
* Bloom filter: http://www.partow.net/programming/bloomfilter/index.html

# Runnig a Test

* First, clone the above libraries, and the U-PSI file. 
* Second, install the libraries and unzip the Feather file. 
* Third, run the following command lines in order:

      cd Feather-implementation-master
    
      g++ -c Rand.cpp -c Hashtable.cpp -c Polynomial.cpp -c Server.cpp -c Client.cpp
    
      g++  -I$home/homeDirectory/include -I/homeDirectory/bloom_filter/bloom_filter.hpp Rand.o Hashtable.o Polynomial.o Server.o Client.o test.cpp  -o test  -lntl -lgmpxx -lgmp -lcryptopp
    
      ./test
    
    
* Note that in the above "homeDirectory" should be replaced with the name of your machine home directory.
