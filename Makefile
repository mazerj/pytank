
DEST=/auto/share/pypeextra

install:
	chmod +x showtank.py; cp showtank.py $(DEST)
	chmod +x tank2hdf5.py; cp tank2hdf5.py $(DEST)
	chmod +x exper2hdf5.sh; cp exper2hdf5.sh $(DEST)/exper2hdf5
	chmod +x getexpers.sh; cp getexpers.sh $(DEST)/getexpers
	cp *.m $(DEST)

clean:
	@/bin/rm -f *.pyc \#*~
