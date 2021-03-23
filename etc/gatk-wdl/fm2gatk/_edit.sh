#/bin/sh

find . -type f -name '*.wdl' -exec grep -H 'docker: dockerImage' {} \;
find . -type f -name '*.wdl' -exec sed -i -e 's/docker: dockerImage/#docker: dockerImage/g' {} \;

#pip3 install biowdl-input-converter
#pip3 install chunked-scatter
