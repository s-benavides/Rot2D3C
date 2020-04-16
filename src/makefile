# Compiler options
# For the Intel F95 compiler use -c -w -O when compiling. In RH9 systems, 
# you also need libc-2.2.93 and link with -lc-2.2.93. For the Absoft F90 
# compiler, use -c -YEXT_SFX=__ -YEXT_NAMES=LCS -O to compile and set the 
# variable UNDER to 0. When UNDER is equal to 1 you can use the variable 
# APPEND to add characters to external functions before compiling.
# When EXTEN is equal to 1 you can use the variable F90EXT 
# to change the extension of the F90/95 files.

COMP     = mpifort -c
LINK     = mpifort 
UNDER    = 0
EXTEN    = 0
APPEND   = _
F90EXT   = f90
LMPI     = -lmpi -lm#-limf
IDIR     = 
LDIR     =
FFTWLDIR = /home/santiago_b/FFTW/fftlib2/lib
LFFTW    = -ldrfftw -ldfftw -L$(FFTWLDIR) 
SHELL    = /bin/bash
MAKE     = /usr/bin/make
PERL     = /usr/bin/perl
ECHO     = echo
BIN      = ../bin

default:
	@$(ECHO) "Usage: make <action>"
	@$(ECHO) "Possible actions are: hd2D, mhd2D, hall25, clean, dist"

fftp.o:
	$(COMP) $(IDIR) fftp.$(F90EXT)

fftp2D.o:
	@$(MAKE) fftp.o
	$(COMP) $(IDIR) fftp_mod.$(F90EXT)
	$(COMP) $(IDIR) $(LDIR) $(LMPI) $(LFFTW) fftp2D.$(F90EXT)

pseudospec2D_inc.o:
	$(COMP) $(IDIR) pseudospec2D_mod.$(F90EXT)
	$(COMP) $(IDIR) pseudospec2D_inc.$(F90EXT)

pseudospec2D_TurbPhi.o:
	$(COMP) $(IDIR) pseudospec2D_mod.$(F90EXT)
	$(COMP) $(IDIR) pseudospec2D_TurbPhi.$(F90EXT)

pseudospec25D_inc.o:
	$(COMP) $(IDIR) pseudospec25D_inc.$(F90EXT)


TurbPhi:
	@$(MAKE) modify
	@$(MAKE) fftp2D.o
	@$(MAKE) pseudospec2D_TurbPhi.o
	$(COMP) $(IDIR) TurbPhi.para.$(F90EXT)
	$(LINK) fftp.o fftp2D.o fftp_mod.o pseudospec2D_mod.o \
	  pseudospec2D_TurbPhi.o TurbPhi.para.o $(LDIR) $(LMPI) $(LFFTW) -o $(BIN)/TurbPhi
	@$(MAKE) undo


modify:
	if [ ${UNDER} -eq 1 ]; then \
	$(PERL) -i.bak -pe 's/fftw_f77\(/fftw_f77${APPEND}\(/g' *.f90 ; \
	for item in `cat external`; do \
	  `echo '$(PERL) -i -pe 's/$$item/$${item}${APPEND}/g' *.f90'`; \
	done; fi
	if [ ${EXTEN} -eq 1 ]; then \
        for file in *.f90; do \
          if [ ${UNDER} -eq 0 ]; then \
             cp $$file $${file}.bak; \
          fi; \
          mv $$file $${file%.f90}.$(F90EXT); \
        done; fi

undo:
	rm -f *.f
	if [ -a fftp.f90.bak ]; then \
        for file in *.f90.bak; do \
          mv $$file $${file%.bak}; \
        done; fi

clean:
	rm -f *.o *.d *.s *.int *.inf *.mod work.pc* ifc*
	@$(MAKE) undo

dist:
	rm -f $(BIN)/TurbPhi* 

sclean: 
	rm -f ../outs/spectrum.*.txt *bal.txt misc.txt ../outs/*.out ../ins/*.out ../outstransfer*.txt time_*.txt ../outs/fluxes.*.txt
	@$(MAKE) clean
