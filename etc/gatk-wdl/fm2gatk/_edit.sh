#/bin/sh

find . -type f -exec grep -H 'docker: dockerImage' {} \;
find . -type f -exec sed -i -e 's/docker: dockerImage/#docker: dockerImage/g' {} \;

