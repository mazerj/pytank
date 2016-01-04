
DEST=/auto/share/pypeextra

install:
	chmod +x tank2hdf5.py; cp tank2hdf5.py $(DEST)
	chmod +x exper2pyt.sh; cp exper2pyt.sh $(DEST)/exper2pyt
	chmod +x getexpers.sh; cp getexpers.sh $(DEST)/getexpers
	cp *.m $(DEST)

clean:
	@/bin/rm -f *.pyc \#*~
