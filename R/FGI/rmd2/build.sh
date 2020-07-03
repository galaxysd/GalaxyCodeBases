#!/usr/bin/env sh
#
# For macOS, `ln -s /usr/local/lib/R/site-library/littler/bin/r /usr/local/bin/littler`
# For Linux, `sudo ln -s /usr/lib/R/site-library/littler/bin/r /usr/local/bin/littler`

#rmarkdown::render('20200517rep.Rmd',output_format='all',output_dir='./out')
#rmarkdown::render('20200608rep.Rmd',output_format='all',output_dir='./out')
Rscript --vanilla -e "rmarkdown::render('20200608rep.Rmd',output_format='all',output_dir='./out')"
Rscript --vanilla -e "rmarkdown::render('20200610rep.Rmd',output_format='all',output_dir='./out')"
