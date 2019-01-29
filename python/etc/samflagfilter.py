#!/usr/bin/env python3

def main():
    import sys
    print("This is the name of the script: ",sys.argv[0])
    print("Number of arguments: ", len(sys.argv))
    print("The arguments are: " , str(sys.argv))
    if len(sys.argv) == 1 :
        print('[x]No arguments given!',file=sys.stderr,flush=True)
        exit(0)
    print('Reading',sys.argv[1])

if __name__ == "__main__":
    main()
