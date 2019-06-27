# Instruction

First, install some of the requirements provided as packages, as below.

Then, install *remained* requirements.

# Requirements

 * Linux CentOS 6.9
 * Perl 5.22.0
 * Python 3.4
 * cutadapt 2.3 <https://pypi.org/project/cutadapt/>
 * bwa 0.7.17
 * samtools 1.9
 * bcftools 1.9
 * R 3.6.0
 * Son of Grid Engine 8.1.9 <https://arc.liv.ac.uk/downloads/SGE/releases/8.1.9/>

# Package Install

## CentOS 6

```bash
sudo yum install https://dl.fedoraproject.org/pub/epel/epel-release-latest-6.noarch.rpm
sudo yum install centos-release-scl
sudo yum install rh-python34
sudo yum install rh-perl524 bwa perl
sudo yum --disablerepo='epel' groupinstall 'Development Tools'
sudo yum install https://arc.liv.ac.uk/downloads/SGE/releases/8.1.9/gridengine-8.1.9-1.el6.x86_64.rpm

pip3 install cutadapt
```

See <https://arc.liv.ac.uk/downloads/SGE/releases/8.1.9/> for details of GridEngine install.

## CentOS 7

```bash
sudo yum install https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm
sudo yum install centos-release-scl
sudo yum install rh-python35
sudo yum install rh-perl524 bwa perl
sudo yum --disablerepo='epel' groupinstall 'Development Tools'
sudo yum install https://copr-be.cloud.fedoraproject.org/results/loveshack/SGE/epel-7-x86_64/00756477-gridengine/gridengine-8.1.9-2.el7.x86_64.rpm

pip3 install cutadapt
```

See <https://copr.fedorainfracloud.org/coprs/loveshack/SGE/> for details of GridEngine install.

## Debian 9 (Stretch)

Follow <https://backports.debian.org/Instructions/> to add `deb http://deb.debian.org/debian stretch-backports main` to `/etc/apt/sources.list.d/`.

```bash
sudo apt install perl python3 r-base build-essential gridengine-client
sudo apt install -t stretch-backports bwa samtools bcftools

pip3 install cutadapt
```
