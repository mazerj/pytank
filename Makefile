
DEST=/auto/share/pypeextra

# set HDF5DUMP=. to dump/load from current dir
HDF5DUMP=/auto/th5

install:
	chmod +x tank2hdf5.py; cp tank2hdf5.py $(DEST)
	chmod +x exper2pyt.sh; cp exper2pyt.sh $(DEST)/exper2pyt
	chmod +x getexpers.sh; cp getexpers.sh $(DEST)/getexpers
	cp pyt_*.m $(DEST)
	cp stab_*.m $(DEST)
	cp list_tdtblocks $(DEST)
	cp tdtsnipcolors.m $(DEST)
	sed -i "s^HDF5DUMP^\'$(HDF5DUMP)\'^" $(DEST)/tank2hdf5.py
	sed -i "s^HDF5DUMP^\'$(HDF5DUMP)\'^" $(DEST)/pyt_load.m

test:

clean:
	@/bin/rm -f *.pyc \#*~
