g++ soap.coverage.cpp -o soap.coverage -lz -pthread -lboost_regex -O2 -fomit-frame-pointer -march=native -funroll-loops
#g++ soap.coverage.cpp -o soap.coverage.static -lz -pthread -l"$regexlib" -O2 -fomit-frame-pointer -march=native -funroll-loops -static
