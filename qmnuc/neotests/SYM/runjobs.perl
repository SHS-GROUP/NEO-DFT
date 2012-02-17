#!/usr/bin/perl

  #############################################################################
  #                                                                           #
  #                           --- runjobs.perl ---                            #
  #                                                                           #
  #                   A script to run GAMESS jobs in batch                    #
  #                                                                           #
  #                                                                           #
  #                            August 2003 by CWS                             #
  #                                                                           #
  #############################################################################


#
#   put all jobs in a hash
#

    %files = ( 1, "sym01",
               2, "sym02",
               3, "sym03",
               4, "sym04",
               5, "sym05",
               6, "sym06",
               7, "sym07",
               8, "sym08",
               9, "sym09" );

    %final_energies = ( 1, "-75.9377419108",
                        2, "-75.9377419106",
			3, "-225.3407001554",
			4, "-55.3317117774",
			5, "-39.5571307715",
			6, "-265.5106745169",
			7, "-1.0481625100",
			8, "-1.0481625096",
		        9, "-919.5868399758" );

#
#   now run them, collect the energies, and print to a log file
#

    $search_string = "FINAL RHF ENERGY IS";
    $search_string2 = "exited gracefully";

    @energy = ();
    
    open( DATA, ">runjobs".".log" );
	 
    for( $i = 1 ; $i <= 9 ; ++$i ) {

         $job     = $files{$i};
         $log     = "$files{$i}".".log";
         $command = "/home/simon/gamess/rungms $job >& $log ";

	 if( ! -e $log ) { system( $command ); }  

         open( FILETOREAD, "<$files{$i}".".log" );

         while( $line = <FILETOREAD> ) {
         
		 if( $line =~ /$search_string/ ) { 
                     
                     @line_chunks = split( /\s+/, $line );
		     $energy[$i]  = $line_chunks[5]; } 

                 if( $line =~ /$search_string2/ ) {
                                                                                
                     print "job $i has ended gracefully\n";
                     $j = $j + 1; } }

         print DATA "$final_energies{$i}\n";
	 print DATA "$energy[$i] \t $files{$i}.log\n\n"; }

#--------- warn the user if certain jobs did not end gracefully
                                                                                
                 if( $j != 9 ) {
                                                                                
            print "one or more of the test jobs did not end properly\n"; }

    exit
