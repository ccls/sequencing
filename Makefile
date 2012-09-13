
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


#	apparently by tagging these make targets as ".PHONY"
#	they will always run when called and not hide behind ...
#	make: Nothing to be done for `all'.
#		OR
#	make: `app1' is up to date.
#	If files are actually created, this is probably not
#	necessary as there are file dates to compare.
#.PHONY: all blat blast bowtie bowtie2
#.PHONY: blat blast bowtie bowtie2 trinity

#	all is the default target that is used when none are given
all: blat bowtie bowtie2 blast trinity

#install: all
install:
	@echo "INSTALLING ALL"
	$(MKDIR) $(BASE_BIN_DIR)
	cd $(BOWTIE) && cp bowtie bowtie-build bowtie-inspect $(BASE_BIN_DIR)
	cd $(BOWTIE2) && cp bowtie2 bowtie2-build bowtie2-align bowtie2-inspect $(BASE_BIN_DIR)
	cp $(BLAT)/bin/* $(BASE_BIN_DIR)
	cd $(BLAST) && make install
	cp -r $(TRINITY)/* $(BASE_BIN_DIR)
	cd rins_core && cp *pl $(BASE_BIN_DIR)
#	mv $(BASE_BIN_DIR)/rins.pl $(BASE_BIN_DIR)/rins.pl.original
#	mv $(BASE_BIN_DIR)/blastn_cleanup.pl $(BASE_BIN_DIR)/blastn_cleanup.pl.original
#	mv $(BASE_BIN_DIR)/write_result.pl $(BASE_BIN_DIR)/write_result.pl.original
	cp rins.pl blastn_cleanup.pl write_result.pl $(BASE_BIN_DIR)
#	mv $(BASE_BIN_DIR)/$(TRINITY)/util/SAM_filter_out_unmapped_reads.pl $(BASE_BIN_DIR)/$(TRINITY)/util/SAM_filter_out_unmapped_reads.pl.original
	cp SAM_filter_out_unmapped_reads.pl $(BASE_BIN_DIR)/$(TRINITY)/util/
#	mv $(TRINITY)/util/SAM_filter_out_unmapped_reads.pl $(TRINITY)/util/SAM_filter_out_unmapped_reads.pl.original
#	cp SAM_filter_out_unmapped_reads.pl $(TRINITY)/util/
	@echo "DONE INSTALLING ALL"
	@echo
	@echo "Add  $(BASE_BIN_DIR) TO YOUR PATH"
	@echo

blat:
	@echo "MAKING BLAT"
	cd $(BLAT) && $(MKDIR) bin && make C_INCLUDE_PATH=/opt/local/include/ PNGLIB=/opt/local/lib/libpng.a BINDIR=`pwd`/bin

bowtie:
	@echo "MAKING BOWTIE"
	cd $(BOWTIE) && make

bowtie2:
	@echo "MAKING BOWTIE2"
	cd $(BOWTIE2) && make

blast:
	@echo "MAKING BLAST"
	#	Use BASE_DIR as blast creates bin/, lib/ and include/
	cd $(BLAST) && ./configure --prefix=$(BASE_DIR) && make

trinity:
	@echo "MAKING TRINITY"
	cd $(TRINITY) && make

#rins:
#	@echo "MAKING RINS WOULD BE NICE, BUT JUST SCRIPTS"
#
#ccls:
#	@echo "MAKING CCLS WOULD BE NICE, BUT JUST SCRIPTS"



clean:
	@echo "sparkling"
	cd $(BLAT) && make clean && rm -f */*/*.a
	/bin/rm -rf $(BLAST)/*-Debug*
	cd $(BOWTIE) && make clean
	cd $(BOWTIE2) && make clean
	cd $(TRINITY) && make clean
#	rins nothing to clean
#	ccls nothing to clean

test:
	@echo "testing is nice, but not today"
#	cd $(BLAT) && make test
#	cd $(BLAST) && make test
#	cd $(BOWTIE) && make test
#	cd $(BOWTIE2) && make test
#	cd $(TRINITY) && make test
#	rins
#	ccls

