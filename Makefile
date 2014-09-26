
DEST=/auto/share/pypeextra

install:
	chmod +x showtank.py; cp showtank.py $(DEST)
	chmod +x tank2hdf5.py; cp tank2hdf5.py $(DEST)
	chmod +x exper2hdf5; cp exper2hdf5 $(DEST)
	cp *.m $(DEST)



clean:
	/bin/rm *.pyc \#*~
