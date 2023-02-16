
==============================
  T0 -> 0.0-2.63
  T1 -> 0-1
  T2 -> 2
------------------------------
 -0- /home/export/online1/mdt00/shisuan/swyaow/luhao/lh/speed_test/Vina@QNLM/build/linux/release/./vinardo 
   -1- slave__Waiting_For_Task ([192/192] T0 0.63/v59535 0.62/v59535 0.61/v59535...)
   -1- main ([0/3])
     -2- PMPI_Recv ([2/2] T1 0/v59535 1/v59535)
     -2- main_ori ([0/1])
       -3- std::vector<std::__cxx11::basic_string<char, ([0/1])
         -4- operator ([0/1])
           -5- free ([0/1])
             -6- unlocked_free ([0/1])
               -7- signal handler called ([0/1])
                 -8- swch_catchsig ([0/1])
                   -9- print_stack ([0/1])
                     -10- __backtrace ([0/1])
                       -11- _Unwind_Backtrace ([0/1])
                         -12- uw_init_context_1 ([0/1])
                           -13- uw_frame_state_for ([0/1])
                             -14- _Unwind_Find_FDE  at /usr1/comp9/swgcc710/gcc-7.1.0/libgcc/unwind-dw2-fde-dip.c:458 ([0/1])
                               -15- _Unwind_Find_registered_FDE ([0/1])
                                 -16- search_object  at /usr1/comp9/swgcc710/gcc-7.1.0/libgcc/unwind-dw2-fde.c:989 ([0/1])
                                   -17- init_object  at /usr1/comp9/swgcc710/gcc-7.1.0/libgcc/unwind-dw2-fde.c:799 ([0/1])
                                     -18- start_fde_sort ([0/1])
                                       -19- malloc ([0/1])
                                         -20- __pthread_mutex_lock ([0/1])
                                           -21- __lll_lock_wait  at lowlevellock.c:45 ([1/1] T2 2/v59535)
==============================

==========================================================
node(taskid):svrstart,wait,sout
vn059535(0       ):	0.45	10.46	10.46
before:0.006435,scan:31.009445,first_data:30.009160,process:1.000427,show:0.020656,total:31.036678

