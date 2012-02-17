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

    %files = ( 1, "diropn01",
               2, "diropn02",
               3, "diropn03",
               4, "diropn04",
               5, "diropn05",
               6, "diropn06",
               7, "diropn07",
               8, "diropn08",
               9, "diropn09",
              10, "diropn10",
              11, "diropn11",
              12, "diropn12",
              13, "diropn13",
              14, "diropn14",
              15, "diropn15" );

    %final_energies = ( 1, "-75.4589899217",
                        2, "-75.4589899216",
			3, "-75.4589899218",
			4, "-75.4589899216",
			5, "-75.4589899216",
			6, "-75.4884497844",
			7, "-75.4884497844",
			8, "-75.4884497844",
			9, "-75.4884497844",
		       10, "-75.4884497844",
		       11, "-75.4884497844",
		       12, "-75.4884497844",
		       13, "-75.4884497844",
		       14, "-75.4884497844",
		       15, "-75.4884497844" );
		       
#
#   now run them, collect the energies, and print to a log file
#

    $search_string  = "FINAL ROHF ENERGY IS";
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
