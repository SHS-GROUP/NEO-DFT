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

    %files = ( 1, "dirc01",
               2, "dirc02",
               3, "dirc03",
               4, "dirc04",
               5, "dirc05",
               6, "dirc06",
               7, "dirc07",
               8, "dirc08",
               9, "dirc09",
              10, "dirc10",
              11, "dirc11",
              12, "dirc12",
              13, "dirc13",
              14, "dirc14",
              15, "dirc15" );

    %final_energies = ( 1, "-515.8466736733",
                        2, "-515.8466736819",
			3, "-515.8466736818",
			4, "-515.8466736819",
			5, "-515.8466736823",
			6, "-515.9249079568",
			7, "-515.9249079653",
			8, "-515.9249079656",
			9, "-515.9249079653",
		       10, "-515.9249079654",
		       11, "-515.9249079568",
		       12, "-515.9249079656",
		       13, "-515.9249079656",
		       14, "-515.9249079653",
		       15, "-515.9249079654" );

#
#   now run them, collect the energies, and print to a log file
#

    $search_string  = "FINAL RHF ENERGY IS";
    $search_string2 = "exited gracefully";

    @energy = ();
    
    open( DATA, ">runjobs".".log" );
	 
    for( $i = 1 ; $i <= 15 ; ++$i ) {

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

                 if( $j != 15 ) { 

            print "one or more of the test jobs did not end properly\n"; }
    exit

