
BASE_DIR     = ${HOME}/RINS_BASE
BASE_BIN_DIR = ${BASE_DIR}/bin
BLAT    = blatSrc34
BLAST   = ncbi-blast-2.2.27+-src/c++
BOWTIE  = bowtie-0.12.8
BOWTIE2 = bowtie2-2.0.0-beta7
TRINITY = trinityrnaseq_r2012-06-08

#	mkdir will raise error if dir exists
#	mkdir -p will not and would create full path
MKDIR        = mkdir -p


#	the @ prefix means the line will be executed, but not printed


#	Do I need the &&s?  Why not just multiple lines?
#	Each line is its own thing so if a cd is used, needs to be same line.


#	all is the default target that is used when none are given
all: blat bowtie bowtie2 blast trinity

#install: all
install: install_all install_bowtie install_bowtie2 install_blat install_blast install_trinity install_scripts
	@printf "\nDONE INSTALLING ALL\n\n"
	@printf "Add  $(BASE_BIN_DIR) TO YOUR PATH\n\n"

install_scripts:
	@printf "\nINSTALLING SCRIPTS\n\n"
	cp rins_core/*pl $(BASE_BIN_DIR)
	cp rins_ccls/*pl $(BASE_BIN_DIR)
	cp rins_ccls/*rb $(BASE_BIN_DIR)
	cp SAM_filter_out_unmapped_reads.pl $(BASE_BIN_DIR)/util/

install_all:
	@printf "\nINSTALLING ALL\n\n"
	$(MKDIR) $(BASE_BIN_DIR)


blat:
	@printf "\nMAKING BLAT\n\n"
	cd $(BLAT) && $(MKDIR) bin && make C_INCLUDE_PATH=/opt/local/include/ PNGLIB=/opt/local/lib/libpng.a BINDIR=`pwd`/bin

install_blat:
	@printf "\nINSTALLING BLAT\n\n"
	cp $(BLAT)/bin/* $(BASE_BIN_DIR)

clean_blat:
	@printf "\nCLEANING BLAT\n\n"
	cd $(BLAT) && make clean && rm -f */*/*.a


bowtie:
	@printf "\nMAKING BOWTIE\n\n"
	cd $(BOWTIE) && make

install_bowtie:
	@printf "\nINSTALLING BOWTIE\n\n"
	cd $(BOWTIE) && cp bowtie bowtie-build bowtie-inspect $(BASE_BIN_DIR)

bowtie2:
	@printf "\nMAKING BOWTIE2\n\n"
	cd $(BOWTIE2) && make

install_bowtie2:
	@printf "\nINSTALLING BOWTIE2\n\n"
	cd $(BOWTIE2) && cp bowtie2 bowtie2-build bowtie2-align bowtie2-inspect $(BASE_BIN_DIR)

blast:
	@printf "\nMAKING BLAST\n\n"
#	Use BASE_DIR as blast creates bin/, lib/ and include/
	cd $(BLAST) && ./configure --prefix=$(BASE_DIR) && make

install_blast:
	@printf "\nINSTALLING BLAST\n\n"
	cd $(BLAST) && make install

trinity:
	@printf "\nMAKING TRINITY\n\n"
	cd $(TRINITY) && make

install_trinity:
	@printf "\nINSTALLING TRINITY\n\n"
#	this works, but copies too much
#	cp -r $(TRINITY)/* $(BASE_BIN_DIR)
	cp -r $(TRINITY)/Analysis $(BASE_BIN_DIR)
	cp -r $(TRINITY)/Butterfly $(BASE_BIN_DIR)
	cp -r $(TRINITY)/Chrysalis $(BASE_BIN_DIR)
	cp -r $(TRINITY)/Inchworm $(BASE_BIN_DIR)
	cp -r $(TRINITY)/PerlLib $(BASE_BIN_DIR)
	cp -r $(TRINITY)/PerlLibAdaptors $(BASE_BIN_DIR)
	cp -r $(TRINITY)/Trinity.pl $(BASE_BIN_DIR)
	cp -r $(TRINITY)/trinity-plugins $(BASE_BIN_DIR)
	cp -r $(TRINITY)/util $(BASE_BIN_DIR)

#LICENSE
#Makefile
#README
#Release.Notes
#docs
#misc
#notes


#rins:
#	@printf "\nMAKING RINS WOULD BE NICE, BUT JUST SCRIPTS\n\n"
#
#ccls:
#	@printf "\nMAKING CCLS WOULD BE NICE, BUT JUST SCRIPTS\n\n"



clean: clean_blat
	@printf "\nCLEANING\n\n"
	/bin/rm -rf $(BLAST)/*-Debug*
	cd $(BOWTIE) && make clean
	cd $(BOWTIE2) && make clean
	cd $(TRINITY) && make clean
#	rins nothing to clean
#	ccls nothing to clean

test:
	@printf "\ntesting is nice, but not today\n\n"
#	cd $(BLAT) && make test
#	cd $(BLAST) && make test
#	cd $(BOWTIE) && make test
#	cd $(BOWTIE2) && make test
#	cd $(TRINITY) && make test
#	rins
#	ccls

