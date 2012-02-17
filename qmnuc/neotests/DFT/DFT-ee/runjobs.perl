#!/usr/bin/perl

  #############################################################################
  #                                                                           #
  #                                                                           #
  #############################################################################

#
#   put all jobs in a hash
#
    %files = ( 1, "FHF_B3LYP_DFT-ee_aug-cc-pVTZ",
               2, "FHF_LYP_DFT-ee_aug-cc-pVTZ",
               3, "H2O_C2v_B3LYP_DFT-ee" );
               

    %dft_energies   = ( 1, "-200.3460818153",
                        2, "-200.2802816520",
                        3, "-76.2984879219" );

    %ee_corr        = ( 1, "-17.3325608427",
                        2, "-0.7351066832",
                        3, "-7.4981684235" );



#
#   now run them, collect the energies, and print to a log file
#

    $search_string1 = 'FINAL';
    $search_string2 = 'DFT EXCHANGE';
    $search_string5 = "exited gracefully";

    @dftenergy = ();
    @eecorr = ();
    
    open( DATA, ">runjobs".".log" );
    print DATA  "\n FIRST LINE IS ENERGY COMPONENT FROM TEST JOB LOG FILE \n";
    print DATA " SECOND LINE IS CORRECT VALUE \n \n \n";

    $j=0;	 
    for( $i = 1 ; $i <= 3 ; ++$i ) {

         $job     = $files{$i};
         $log     = "$files{$i}".".log";
         $command = "/home/simon/gamess/rungms $job >& $log ";

	 if( ! -e $log ) { system( $command ); }  

         open( FILETOREAD, "<$files{$i}".".log" );

         while( $line = <FILETOREAD> ) {
         
		 if( $line =~ /$search_string1/ ) {
                     @line_chunks = split( /\s+/, $line );
		     $dftenergy[$i]  = $line_chunks[5]; } 

		 if( $line =~ /$search_string2/ ) { 
                     @line_chunks = split( /\s+/, $line );
		     $eecorr[$i]  = $line_chunks[7]; } 

                 if( $line =~ /$search_string5/ ) {

                     print "job $i has ended gracefully\n";
                     $j = $j + 1; } }

         print DATA "  $files{$i}\n";
	 print DATA " DFT E0 = $dftenergy[$i] \n"; 
         print DATA " DFT E0 = $dft_energies{$i} \n \n"; 
	 print DATA " EXCH + CORR = $eecorr[$i] \n"; 
         print DATA " EXCH + CORR = $ee_corr{$i} \n \n"; 
         print DATA "********************************************** \n \n \n "; }

#--------- warn the user if certain jobs did not end gracefully

                 if( $j != 3 ) { 
            print "one or more of the test jobs did not end properly\n"; }
                 if( $j = 3 ) { 
            print "all of the test jobs exited gracefully\n"; }
    exit

