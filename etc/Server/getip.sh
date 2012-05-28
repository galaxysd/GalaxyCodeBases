#!/bin/bash

# http://www.linuxquestions.org/questions/linux-networking-3/how-can-i-get-my-external-ip-address-from-behind-a-nat-333878/

get-current-ip() {

xc=1
until (( xc == 0 )) ; do

     case $(( RANDOM % 5 )) in

          0) ip=$(wget -t 2 -T 5 -q -O- http://showip.codebrainz.ca/) ;;
#         1) ip=$(wget -t 2 -T 5 -q -O- http://www.whatismyip.com/automation/n09230945.asp) ;;
          2) ip=$(wget -t 2 -T 5 -q -O- http://www.showmyip.com/simple/) ;;
          3) ip=$(wget -t 2 -T 5 -q -O- http://cfaj.freeshell.org/ipaddr.cgi) ;;
          4) ip=$(wget -t 2 -T 5 -q -O- https://secure.informaction.com/ipecho/) ;;
          1) ip=$(wget -t 2 -T 5 -q -O- http://icanhazip.com/) ;;

     esac

     xc=$?

done

#echo -n "${ip//[^0-9.]}"   #if you don't want a trailing newline
echo "$ip"                #if you want or don't mind newlines
}

get-current-ip

