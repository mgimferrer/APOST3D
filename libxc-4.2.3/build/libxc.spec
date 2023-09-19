# RPM spec file for libxc.
# This file is used to build Redhat Package Manager packages for the
# libxc.  Such packages make it easy to install and uninstall
# the library and related files from binaries or source.
#
# This spec file is for version 4.2.3 of libxc; the appropriate
# version numbers are automatically substituted in to libxc.spec.in
# by the configure script.  However, libxc.spec.in may need to be
# modified for future releases, if the list of installed files
# or build commands change.
#
# RPM.  To build, use the command: rpm --clean -ba libxc.spec
#
# Alternatively, you can just use 'make rpm'.
#
Name: libxc
Summary: Library of exchange and correlation functionals to be used in DFT codes
Version: 4.2.3
Release: 1
Provides: %{name}
License: MPL 2.0
Group: Applications/Scientific
Prefix: /usr
BuildRoot: %{_tmppath}/%{name}-%{version}-buildroot
Source: http://www.tddft.org/programs/octopus/download/%{name}-%{version}.tar.gz
URL: http://www.tddft.org/programs/octopus/wiki/index.php/Libxc

%description 
Libxc is a library of exchange and correlation functionals. Its
purpose is to be used in codes that implement density-functional
theory. For the moment, the library includes most of the local density
approximations (LDAs), generalized density approximation (GGAs), and
meta-GGAs. The library provides values for the energy density and its
1st, 2nd, and (for the LDAs) 3rd derivatives.

%prep
rm -rf $RPM_BUILD_ROOT
%setup -q

# The installation is also performed in the %%build stage because the
# code has to be configured twice, with and without MPI support, and
# cleaned in between.
%build
%configure \
  CC="icc" \
  CPP="icc -E" \
  FC="ifort" \
  FCFLAGS="-u -fpp1 -nbs -pc80 -pad -align -unroll-aggressive -O3 -ip -no-fp-port -mno-ieee-fp -vec-report0 -no-prec-div -parallel -qopenmp -qmkl" \
  CFLAGS="-O3" \
  CPPFLAGS="" \
  LDFLAGS="-parallel -qopenmp -qmkl" \
  --disable-static

make

make install DESTDIR=${RPM_BUILD_ROOT}

%clean
rm -rf ${RPM_BUILD_ROOT}


%post


%preun


%files
%defattr(-,root,root,0755)
%doc README NEWS COPYING AUTHORS ChangeLog
%{_libdir}/*
%{_includedir}/*
