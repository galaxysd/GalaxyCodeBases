#/bin/sh

find . -type f -name '*.wdl' -exec grep -H 'docker: dockerImage' {} \;
find . -type f -name '*.wdl' -exec sed -i -e 's/docker: dockerImage/#docker: dockerImage/g' {} \;

