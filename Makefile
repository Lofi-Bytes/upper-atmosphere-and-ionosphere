iono.so: iono.f90 iono.pyf
	f2py -c --fcompiler=gfortran iono.pyf iono.f90

iono.pyf: iono.f90
	f2py -m iono -h iono.pyf $<

clean:
	rm -f iono.so *~ iono.pyf

iono:iono.so
	python iono.py