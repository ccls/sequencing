
BASE_DIR     = ${HOME}/RINS_BASE
BASE_BIN_DIR = ${BASE_DIR}/bin
BLAT    = blatSrc34
BLAST   = ncbi-blast-2.2.27+-src/c++
BOWTIE  = bowtie-0.12.8
BOWTIE2 = bowtie2-2.0.0-beta7
#TRINITY = trinityrnaseq_r2012-06-08
TRINITY = trinityrnaseq_r2012-10-05
BWA = bwa-0.6.2
PRICE = PriceSource120527
RAY = Ray-20121031
VELVET = Velvet-20121031

#	mkdir will raise error if dir exists
#	mkdir -p will not and would create full path
MKDIR        = mkdir -p


#	In the preferred compilation order ...
TARGETS = blat bowtie bowtie2 blast bwa trinity price velvet ray


#	the @ prefix means the line will be executed, but not printed


#	Do I need the &&s?  Why not just multiple lines?
#	Each line is its own thing so if a cd is used, needs to be same line.


#	all is the default target that is used when none are given
#all: blat bowtie bowtie2 blast bwa trinity price
all: make-all $(TARGETS)
	@printf "\nDONE MAKING ALL\n\n"

make-all:
	@printf "\nMAKING ALL\n\n"

#install: all
#install: install_all install_bowtie install_bowtie2 install_blat install_blast install_trinity install_bwa install_price install_scripts
install: install-all $(TARGETS:%=install-%) install-scripts
	@printf "\nDONE INSTALLING ALL\n\n"
	@printf "Add  $(BASE_BIN_DIR) TO YOUR PATH\n\n"

install-all:
	@printf "\nINSTALLING ALL\n\n"
	$(MKDIR) $(BASE_BIN_DIR)

install-scripts:
	@printf "\nINSTALLING SCRIPTS\n\n"
	cp rins_core/*pl $(BASE_BIN_DIR)
#	no more perl scripts here 
#	Oops.  Had to modify compress.pl
	cp rins_ccls/*pl $(BASE_BIN_DIR)
	cp rins_ccls/*rb $(BASE_BIN_DIR)
#	using latest trinity which may not need this
#	cp SAM_filter_out_unmapped_reads.pl $(BASE_BIN_DIR)/util/


#clean: clean_blat clean_bwa clean_trinity clean_price
clean: clean-all $(TARGETS:%=clean-%)
	@printf "\nDONE CLEANING\n\n"

clean-all:
	@printf "\nCLEANING ALL\n\n"
#	rins nothing to clean
#	ccls nothing to clean


blat:
	@printf "\nMAKING BLAT\n\n"
	cd $(BLAT) && $(MKDIR) bin && make C_INCLUDE_PATH=/opt/local/include/ PNGLIB=/opt/local/lib/libpng.a BINDIR=`pwd`/bin

install-blat:
	@printf "\nINSTALLING BLAT\n\n"
	cp $(BLAT)/bin/* $(BASE_BIN_DIR)

clean-blat:
	@printf "\nCLEANING BLAT\n\n"
	cd $(BLAT) && make clean && rm -f */*/*.a



bowtie:
	@printf "\nMAKING BOWTIE\n\n"
	cd $(BOWTIE) && make

install-bowtie:
	@printf "\nINSTALLING BOWTIE\n\n"
	cd $(BOWTIE) && cp bowtie bowtie-build bowtie-inspect $(BASE_BIN_DIR)

clean-bowtie:
	@printf "\nCLEANING BOWTIE\n\n"
	cd $(BOWTIE) && make clean



bowtie2:
	@printf "\nMAKING BOWTIE2\n\n"
	cd $(BOWTIE2) && make

install-bowtie2:
	@printf "\nINSTALLING BOWTIE2\n\n"
	cd $(BOWTIE2) && cp bowtie2 bowtie2-build bowtie2-align bowtie2-inspect $(BASE_BIN_DIR)

clean-bowtie2:
	@printf "\nCLEANING BOWTIE2\n\n"
	cd $(BOWTIE2) && make clean


blast:
	@printf "\nMAKING BLAST\n\n"
#	Use BASE_DIR as blast creates bin/, lib/ and include/
	cd $(BLAST) && ./configure --prefix=$(BASE_DIR) && make

install-blast:
	@printf "\nINSTALLING BLAST\n\n"
	cd $(BLAST) && make install

clean-blast:
	@printf "\nCLEANING BLAST\n\n"
	/bin/rm -rf $(BLAST)/*-Debug*



trinity:
	@printf "\nMAKING TRINITY\n\n"
	cd $(TRINITY) && make

install-trinity:
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

clean-trinity:
	@printf "\nCLEANING TRINITY\n\n"
	cd $(TRINITY) && make clean

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

bwa:
	@printf "\nMAKING BWA\n\n"
	cd $(BWA) && make

install-bwa:
	@printf "\nINSTALLING BWA\n\n"
	cp $(BWA)/bwa $(BASE_BIN_DIR)

clean-bwa:
	@printf "\nCLEANING BWA\n\n"
	cd $(BWA) && make clean



price:
	@printf "\nMAKING PRICE\n\n"
	cd $(PRICE) && make

install-price:
	@printf "\nINSTALLING PRICE\n\n"
	cp $(PRICE)/PriceTI $(BASE_BIN_DIR)

clean-price:
	@printf "\nCLEANING PRICE\n\n"
	cd $(PRICE) && make clean



#	The dir name is ray so 'make ray' will not run.
.PHONY: ray
ray:
	@printf "\nMAKING RAY\n\n"
	cd $(RAY) && make

install-ray:
	@printf "\nINSTALLING RAY\n\n"
	cp $(RAY)/Ray $(BASE_BIN_DIR)

clean-ray:
	@printf "\nCLEANING RAY\n\n"
	cd $(RAY) && make clean




velvet:
	@printf "\nMAKING VELVET\n\n"
	cd $(VELVET) && make

install-velvet:
	@printf "\nINSTALLING VELVET\n\n"
	cp $(VELVET)/velvet? $(BASE_BIN_DIR)

clean-velvet:
	@printf "\nCLEANING VELVET\n\n"
	cd $(VELVET) && make clean







test:
	@printf "\ntesting is nice, but not today\n\n"
#	cd $(BLAT) && make test
#	cd $(BLAST) && make test
#	cd $(BOWTIE) && make test
#	cd $(BOWTIE2) && make test
#	cd $(TRINITY) && make test
#	rins
#	ccls

