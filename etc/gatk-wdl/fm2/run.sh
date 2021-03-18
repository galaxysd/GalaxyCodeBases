#!/bin/sh

echo Begin: `date`

JAVA_OPTS="-Dbackend.providers.Local.config.concurrent-job-limit=16" cromwell run -o cromwell_options_no_user.json -i fm2.json fm2.wdl >fm2.log

echo End: `date`
