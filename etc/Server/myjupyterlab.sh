#!/bin/sh

#micromamba create -n nb
#micromamba install -n nb -c conda-forge jupyterlab ipywidgets python[version='>=3.12'] git[version='>=2.44'] conda   #nodejs "libuv[version='>=1.48']"
######micromamba run -n nb pip3 wheel --no-input --use-pep517 --no-deps -w ./pywhls git+https://github.com/Anaconda-Platform/nb_conda_kernels.git
#micromamba run -n nb pip3 install --no-input --use-pep517 git+https://github.com/Anaconda-Platform/nb_conda_kernels.git

JUPYTERCFGPATH=`micromamba run -n nb jupyter --config-dir`
PASSWORD='Salus107'
HASHEDPW=`printf "from jupyter_server.auth import passwd\nprint(passwd('${PASSWORD}'))"|micromamba run -n nb python3`

printf "
c = get_config()  #noqa
c.ServerApp.allow_remote_access = True
c.ServerApp.ip = '*'
c.ServerApp.open_browser = False
c.ServerApp.port = 10000
c.PasswordIdentityProvider.hashed_password = '${HASHEDPW}'
" > ${JUPYTERCFGPATH}/jupyter_server_config.py

printf '{
  "CondaKernelSpecManager": {
    "kernelspec_path": "--user"
  }
}
' > ${JUPYTERCFGPATH}/jupyter_server_config.json

printf '#!/bin/sh\nPYDEVD_DISABLE_FILE_VALIDATION=1 micromamba run -n nb python3 -Xfrozen_modules=off -m jupyterlab_server\n' > ${MAMBA_ROOT_PREFIX}/envs/nb/bin/nb
chmod +x ${MAMBA_ROOT_PREFIX}/envs/nb/bin/nb
ls -l ${MAMBA_ROOT_PREFIX}/envs/nb/bin/nb

exit

micromamba install -c conda-forge ipykernel ipywidgets r-irkernel -n salus
micromamba run -n nb python -m nb_conda_kernels list

micromamba run -n salus jupyter kernelspec list
micromamba run -n salus python -m ipykernel install --user --name salus --display-name "Python (salus)"
micromamba run -n salusg python -m ipykernel install --user --name salusg --display-name "Python (salusG)"
# IRkernel::installspec(name = 'atacR', displayname = 'atacR')
Rscript -e "IRkernel::installspec(name = 'atacR2', displayname = 'atacR2')"

micromamba run -n nb jupyter kernelspec list

#$ cat ~/micromamba/envs/nb/bin/nb
#!/bin/sh
PYDEVD_DISABLE_FILE_VALIDATION=1 micromamba run -n nb python3 -Xfrozen_modules=off -m jupyterlab_server
