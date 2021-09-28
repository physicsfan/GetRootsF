# FMKMF_F90 set to gfortran
# FMKMF_SFTAG not set: using f90
# FMKMF_SPATH not set: using . 
# FMKMF_LINKOPTS not set: using no link options 
# Main program is froots.f90 
# process_fsource called with arg froots.f90 
# froots.f90 Uses Module test_functions
# froots.f90 Uses Module roots
# Full list of modules in froots.f90: test_functions roots 
# Uses test_functions which is in ./test_functions.f90
# process_fsource called with arg ./test_functions.f90 
# Full list of modules in ./test_functions.f90:  
# Uses roots which is in ./roots.f90
# process_fsource called with arg ./roots.f90 
# Full list of modules in ./roots.f90:  

# ------------------Macro-Defs---------------------
F90=gfortran 

# -------------------End-macro-Defs---------------------------

# Here is the link step 
troots:test_functions.o roots.o test_roots.o 
	 $(F90) -o troots test_functions.o roots.o test_roots.o   

# Here are the compile steps
 
test_functions.o:./test_functions.f90  
	 $(F90) -c ./test_functions.f90 

roots.o:./roots.f90  
	 $(F90) -c ./roots.f90 

test_roots.o:test_roots.f90 test_functions.o roots.o 
	 $(F90) -c test_roots.f90 
# This entry allows you to type " make clean " to get rid of
# all object and module files 
clean:
	rm -f -r f_{files,modd}* *.o *.mod *.M *.d V*.inc *.vo \
	V*.f *.dbg album F.err
  
