## Server

```bash
cd /share/dat0/cromwell
#cp -a ../test/bin .
#ln -s ../fq .
cromwell server
```

## Client

```bash
zip -9r fm2wdl.zip fm2.wdl  tasks/
#cromwell submit --workflow-root /share/dat0/test/ -i fmtest.json /share/dat0/test/fm2.wdl -o cromwell_options_no_user.json -p fm2wdl.zip
cromwell submit -i fmtest.json fm2.wdl -o cromwell_options_no_user.json -p fm2wdl.zip
```
