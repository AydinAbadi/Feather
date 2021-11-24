# Update-Simulation

Update Simulation of the protocol proposed in: https://ieeexplore.ieee.org/abstract/document/7934388/

# Dependencies

* GMP: https://gmplib.org/
* Cryptopp: https://www.cryptopp.com
* NTL: https://www.shoup.net/ntl

# Runnig a Test

* First, clone the above libraries, and the Update-Simulation-master file. Then, install the libraries and unzip the Update-Simulation file. After that:

      cd Update-Simulation-code-master
    
      g++  -c  Rand.cpp -c Hashtable.cpp -c Polynomial.cpp 
    
      g++  -I$home/homeDirectory/include Rand.o Hashtable.o Polynomial.o EO-PSI-Update.cpp  -o test  -lntl -lgmpxx -lgmp -lcryptopp
    
      ./test
    
* In the above, "homeDirectory" should be replaced with the name of your machine home directory.
