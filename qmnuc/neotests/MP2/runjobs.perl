#!/usr/bin/perl

  #############################################################################
  #                                                                           #
  #                                                                           #
  #############################################################################

#
#   put all jobs in a hash
#

    %files = ( 1, "mp2-01",
               2, "mp2-02",
               3, "mp2-03",
               4, "mp2-04",
               5, "mp2-05",
               6, "mp2-06" );

    %hf_energies    = ( 1, "-199.5456898232",
                        2, "-5.7375795891",
			3, "-5.7666820138",
			4, "-199.5456898257",
			5, "-200.0495759763", 
                        6, "-5.7192308547" );

    %ee_energies    = ( 1, "-0.5796098727",
                        2, "-0.0667937890",
			3, "-0.0673089824",
			4, "-0.5796098724",
			5, "-0.5617879883", 
                        6, "-0.0298566574" );

    %ep_energies    = ( 1, "-0.0096428665",
                        2, "-0.0270918080",
			3, "-0.0099134293",
			4, "-0.0096428664",
			5, "-0.0106755498",
                        6, "-0.0208198133" );

    %total_energies = ( 1, "-200.1349425624",
                        2, "-5.8314651861",
			3, "-5.8439044255",
			4, "-200.1349425645",
			5, "-200.6220395145",
                        6, "-5.7699073254" );

#
#   now run them, collect the energies, and print to a log file
#

    $search_string1 = 'E\(0\)=';
    $search_string2 = 'E\(2EE\)=';
    $search_string3 = 'E\(2EP\)=';
    $search_string4 = 'E\(MP2\)=';
    $search_string5 = "exited gracefully";

    @rhfenergy = ();
    @eeenergy = ();
    @epenergy = ();
    @tenergy = ();
    
    open( DATA, ">runjobs".".log" );
    print DATA "FIRST LINE IS ENERGY COMPONENT FROM TEST JOB LOG FILE \n";
    print DATA "SECOND LINE IS CORRECT VALUE \n";

    $j=0;	 
    for( $i = 1 ; $i <= 6 ; ++$i ) {

         $job     = $files{$i};
         $log     = "$files{$i}".".log";
         $command = "/home/simon/gamess/rungms $job >& $log ";

	 if( ! -e $log ) { system( $command ); }  

         open( FILETOREAD, "<$files{$i}".".log" );

         while( $line = <FILETOREAD> ) {
         
		 if( $line =~ /$search_string1/ ) {
                     @line_chunks = split( /\s+/, $line );
		     $rhfenergy[$i]  = $line_chunks[2]; } 

		 if( $line =~ /$search_string2/ ) { 
                     @line_chunks = split( /\s+/, $line );
		     $eeenergy[$i]  = $line_chunks[2]; } 

		 if( $line =~ /$search_string3/ ) { 
                     @line_chunks = split( /\s+/, $line );
		     $epenergy[$i]  = $line_chunks[2]; } 

		 if( $line =~ /$search_string4/ ) { 
                     @line_chunks = split( /\s+/, $line );
		     $tenergy[$i]  = $line_chunks[2]; } 

                 if( $line =~ /$search_string5/ ) {

                     print "job $i has ended gracefully\n";
                     $j = $j + 1; } }

         print DATA "  $files{$i}\n";
	 print DATA "  E0 = $rhfenergy[$i] \n "; 
         print DATA " E0 = $hf_energies{$i} \n \n"; 
	 print DATA "  EE = $eeenergy[$i] \n "; 
         print DATA " EE = $ee_energies{$i} \n \n"; 
	 print DATA "  EP = $epenergy[$i] \n "; 
         print DATA " EP = $ep_energies{$i} \n \n"; 
	 print DATA "  ET = $tenergy[$i] \n  ";
         print DATA "ET = $total_energies{$i} \n \n";
         print DATA "********************************************** \n \n \n "; }

#--------- warn the user if certain jobs did not end gracefully

                 if( $j != 6 ) { 
            print "one or more of the test jobs did not end properly\n"; }
                 if( $j = 6 ) { 
            print "all of the test jobs exited gracefully\n"; }
    exit

