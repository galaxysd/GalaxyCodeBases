%define version 2022.07.07

Summary: Clash Premium RPM Package
Name: clash-premium
Version: %{version}
Release: 2
License: non-free
URL: https://github.com/Dreamacro/clash/releases/tag/premium
Group: System
Packager: Yuuki Galaxy
Source1: clash-linux-amd64-v3-%{version}.gz
Source2: clash.service

%description
A rule-based tunnel in Go.
Enable with `sudo systemctl --system enable clash.service`.

%prep
rm -rf $RPM_BUILD_ROOT
mkdir -p $RPM_BUILD_ROOT/%{_bindir}
cp %{SOURCE1} $RPM_BUILD_ROOT/%{_bindir}/clash.gz
gzip -d $RPM_BUILD_ROOT/%{_bindir}/clash.gz
mkdir -p $RPM_BUILD_ROOT/%{_prefix}/lib/systemd/system/
cp %{SOURCE2} $RPM_BUILD_ROOT/%{_prefix}/lib/systemd/system/
exit

%files
%attr(0755, root, root) /usr/bin/clash
%attr(0644, root, root) /usr/lib/systemd/system/clash.service

