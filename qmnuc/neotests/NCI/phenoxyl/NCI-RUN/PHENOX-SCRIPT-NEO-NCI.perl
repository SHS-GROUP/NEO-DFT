#!/usr/bin/perl

#
#  J. H. Skone  Aug 2007
#
#  Script to compute the NEO-NCI energies as a function of basis function
#  center separation for the phenoxyl/phenol self-exchange reaction.
#
#  You must change the $GMS and $VER variables below. $GMS points to your
#  compiled version of GAMESS and $VER is the version number of your 
#  GAMESS executable (usually default version number is 00)
   $GMS = '/scr/jhskone/gamess-jhs-neo-ver2/rungms';
   $VER = '01';
#
#  Inputs:
#          input file -->   A user defined tag for all files, i.e.
#                           if you input 'test' you get files:
#                           test_0.inp,dat,log
#                           and a data file test.dat
#
#          D-A distance->   The donor acceptor distance in angstroms
#
#          # of steps -->   The number of steps to take in reducing the basis
#                           function center distance starting from H Z COORD
#
#          H Z COORD  -->   The initial (positive) proton basis center coordinate
#                           along the Z - principle axis. The coordinate
#                           for the additional (negative) proton basis center will
#                           be generated in the script. This value will be the
#                           largest separation searched as each additional setp
#                           will be for a smaller bfc separation.
#
#
#
#                The script will also look for files 'orbital_file1' and 'orbital_file2'
#                These files should contain the $VEC group generated from the initial 
#                left and right localized Hartree-Fock runs, respectively. Make sure 
#                both of these files contain the $VEC groups and are placed in the directory
#                where this script is run.
#
#
#  Outputs:
#                Two data files are printed entitled "file"-energies.dat and "file"-overlaps.dat
#                where "file" is the user defined tag for all files given as the first argument
#                on the command line. "file"-energies contains the nuclear basis function center
#                separation, the ground state NEO-NCI energy, the first excited state NEO-NCI
#                energy, the Hartree-Fock energy, the unorthogonalized off-diagonal NEO-NCI 
#                matrix element (H12), the off-diagonal overlap matrix element (S12), and the
#                splitting (in cm-1). The output is printed in table form as:
#                
#                nuclear basis center separation, NEO-NCI E0, NEO-NCI E1,  ROHF E, H12, S12,Splitting
#  
#                "file"-overlaps contains the off-diagonal overlap matrix element and the 
#                corresponding  eletronic and nuclear overlap components. The output is 
#                printed in table form as:
#
#                  nuclear basis center separation, Se12, Sp12,  S12, Splitting
#
#                
#  Testing:      to test the script you should use the following input parameters 
#                input file   -->  PHENOX-TEST
#                D-A distance -->  2.40
#                # steps      -->  5
#                H Z COORD    -->  0.285
#                
#                Submit the script on the command line as follows:
#                PHENOX-SCRIPT-NEO-NCI.perl PHENOX-TEST 2.40 5 0.285 &
#                 
#                Compare your results in PHENOX-TEST-energies.dat and PHENOX-TEST-overlaps.dat with 
#                the results in the PHENOX-NEO-NCI-RUN-energies.dat and PHENOX-NEO-NCI-RUN-overlaps.dat 
#                files, respectively.
#                
#                
#
# check for arguments
#

  unless( @ARGV > 3 ) {
	   
  die " USAGE:  <input file> <D_A> <#steps> <H Z COORD>  \n"; }

#
# Read in all input values:
#

  $input          = $ARGV[0];
  $D_A            = $ARGV[1];
  $num_step       = $ARGV[2];

#
# set the initial basis function center coordinates
#

  $zh        = $ARGV[3];
  $zh2       = -$ARGV[3];

# The basis function center separation is presently reduced by 0.01 angstroms
# at each step. This can be changed by changing the values of $disph1 and 
# $disph2 below.

$disph1 = -0.005;
$disph2 =  0.005;


#
# Read in the electronic molecular orbitals
#
     %RowsFromOrbitalFile = ();
     $RowsInOF = 1;
   
     open( ORBITALS, "orbital_file1" )	 
     || die "Can't open the left localized electronic MO file!!!!!\n";
	 
     while( $line = <ORBITALS> ) {
		 
            $RowsFromOrbitalFile{$RowsInOF} = $line;
            ++$RowsInOF; }
		 
     --$RowsInOF;

     %RowsFromOrbitalFile2 = ();
     $RowsInOF2 = 1;
   
     open( ORBITALS2, "orbital_file2" )	 
     || die "Can't open the right localized electronic MO file!!!!!\n";
	 
     while( $line = <ORBITALS2> ) {
		 
            $RowsFromOrbitalFile2{$RowsInOF2} = $line;
            ++$RowsInOF2; }
		 
     --$RowsInOF2;

#
# Initialize arrays and assign search string
# 

  @line_chunks  = ();
  @es1    = ();
  @es2    = ();
  @h11    = ();
  @h12    = ();
  @S12    = ();
  @dist   = ();
  @line_chunks2  = ();
  @HDIAG  = ();
  $string1 = "ENERG OUT";
  $string2 = "VIBRONIC SPLITTING ";


#
# Open and write header to energy output file
#

 open( MATH_DATAFILE, ">$input"."-energies".".dat" );

 print MATH_DATAFILE "\n";
 print MATH_DATAFILE "\n";
 print MATH_DATAFILE "\n";
 print MATH_DATAFILE "   PHENOXYL PHENOL SYSTEM\n";
 print MATH_DATAFILE "   BASIS SET: 6-31G/QZSPDZD        O-O distance: $D_A\n";
 print MATH_DATAFILE "   All energies are in atomic units except the splitting\n";
 print MATH_DATAFILE "   H-H                                                             Unortho    Overlap     Vibronic\n"; 
 print MATH_DATAFILE "   dist (angs.)    NEO-NCI ST 1      NEO-NCI ST 2      NEO-RHF       H12        S12       Splitting (cm-1)\n";
 print MATH_DATAFILE "   ----------------------------------------------------------------------------------------------------------------\n";
 print MATH_DATAFILE "\n";

#
#  Open and write header to overlap output file.
#

 open( MATH_DATAFILE2, ">$input"."-overlaps".".dat" );

 print MATH_DATAFILE2 "\n";
 print MATH_DATAFILE2 "\n";
 print MATH_DATAFILE2 "\n";
 print MATH_DATAFILE2 "    PHENOXYL PHENOL SYSTEM\n";
 print MATH_DATAFILE2 "    BASIS SET: 6-31G/QZSPDZD        O-O distance: $D_A\n";
 print MATH_DATAFILE2 "\n";
 print MATH_DATAFILE2 "    H-H (angs.)\t     Se12 \t     Sp12 \t     St12 \t   Splitting (cm-1)\n";
 print MATH_DATAFILE2 "    -------------------------------------------------------------------------------\n";
 print MATH_DATAFILE2 "\n";
 @line_chunks3  = ();
 @Se12    = ();
 @Sp12    = ();
 $string3 = "SE12 AND SP12";


#
# Begin the energy calculations
#

for( $i = 0; $i <= $num_step; ++$i )  {
	if($zh > 0.025) { 

#
# Generate the left localized NEO-HF input file
#
      open( GAMESSINPFILE, ">$input"."_"."LEFT"."_"."$i".".inp" );
	 
      $LLOC_file_to_run     = "$input"."_"."LEFT"."_"."$i";
      $LLOC_out_file        = "$input"."_"."LEFT"."_"."$i".".log"; 
       system " rm ./$LLOC_file_to_run.dat"; 

	print GAMESSINPFILE "! \n";
	print GAMESSINPFILE "! \n";
	print GAMESSINPFILE "! \n";
	print GAMESSINPFILE "! \n";
	print GAMESSINPFILE "! \n";
	print GAMESSINPFILE "! \n";
	print GAMESSINPFILE "! \n";
        print GAMESSINPFILE " \$CONTRL SCFTYP=ROHF RUNTYP=ENERGY MAXIT=200     \n";
        print GAMESSINPFILE "          INTTYP=HONDO ICHARG=0 MULT=2 \$END\     \n";
        print GAMESSINPFILE " \$SYSTEM TIMLIM=3000 MEMORY=60000000  \$END\     \n";
        print GAMESSINPFILE " \$INTGRL SCHWRZ=.F. \$END\                       \n";
        print GAMESSINPFILE " \$NEO    NEOSCF=DIAGZN BASNUC=QZSPDD NUNIQN=1    \n";
	print GAMESSINPFILE "          IUNIQN(1)=25 EXCH=.T. NUMULT=2 LOCORB=1 \n";
        print GAMESSINPFILE "          NTAUXB=1 NAUXNB=1 IAUXNB(1)=26 \$END\   \n";
	print GAMESSINPFILE " \$BASIS  GBASIS=N31 NGAUSS=6 \$END\              \n";
 	print GAMESSINPFILE " \$GUESS  GUESS=MOREAD NORB=150 \$END\            \n";
	print GAMESSINPFILE " \$DATA\                                          \n";
	print GAMESSINPFILE " \n";
        print GAMESSINPFILE " C1                                          \n";
        print GAMESSINPFILE " O 8.0    0.000000     0.000000    -1.200748 \n"; 
        print GAMESSINPFILE " O 8.0    0.000000     0.000000     1.200748 \n"; 
	print GAMESSINPFILE " C 6.0   -1.114666    -0.400558     1.749899 \n"; 
	print GAMESSINPFILE " C 6.0   -1.307180    -0.207036     3.150282 \n"; 
	print GAMESSINPFILE " C 6.0   -2.479880    -0.616856     3.759596 \n"; 
	print GAMESSINPFILE " C 6.0   -3.498257    -1.226627     3.007797 \n"; 
	print GAMESSINPFILE " C 6.0   -3.323431    -1.426760     1.629474 \n"; 
	print GAMESSINPFILE " C 6.0   -2.153726    -1.033298     1.000684 \n"; 
	print GAMESSINPFILE " C 6.0    1.114666    -0.400558    -1.749899 \n"; 
	print GAMESSINPFILE " C 6.0    2.153727    -1.033309    -1.000695 \n"; 
	print GAMESSINPFILE " C 6.0    3.323415    -1.426764    -1.629516 \n"; 
	print GAMESSINPFILE " C 6.0    3.498237    -1.226624    -3.007839 \n"; 
	print GAMESSINPFILE " C 6.0    2.479865    -0.616833    -3.759629 \n"; 
	print GAMESSINPFILE " C 6.0    1.307176    -0.207021    -3.150281 \n"; 
	print GAMESSINPFILE " H 1.0   -0.510170     0.275234     3.707418 \n"; 
	print GAMESSINPFILE " H 1.0   -2.616289    -0.462710     4.826741 \n"; 
	print GAMESSINPFILE " H 1.0   -4.415323    -1.547947     3.492806 \n"; 
	print GAMESSINPFILE " H 1.0   -4.109267    -1.904896     1.050621 \n"; 
	print GAMESSINPFILE " H 1.0   -2.002195    -1.195765    -0.061871 \n"; 
	print GAMESSINPFILE " H 1.0    2.002199    -1.195787     0.061859 \n"; 
	print GAMESSINPFILE " H 1.0    4.109246    -1.904915    -1.050671 \n"; 
	print GAMESSINPFILE " H 1.0    4.415300    -1.547945    -3.492852 \n"; 
	print GAMESSINPFILE " H 1.0    2.616271    -0.462682    -4.826774 \n"; 
	print GAMESSINPFILE " H 1.0    0.510165     0.275257    -3.707408 \n"; 
	print GAMESSINPFILE " H 1.0    0.000000     0.000000    $zh2      \n"; 
	print GAMESSINPFILE " H 1.0    0.000000     0.000000     $zh      \n"; 
	print GAMESSINPFILE " \$END\ \n";
        for( $k = 1; $k <= $RowsInOF; ++$k ) {
                 
                   print GAMESSINPFILE "$RowsFromOrbitalFile{$k}"; }

	 $command_1 = " $GMS $LLOC_file_to_run $VER >& $LLOC_out_file";
	          system( $command_1 );
#
#   Generate right localized NEO-HF input file
#
      open( GAMESSINPFILE2, ">$input"."_"."RGHT"."_"."$i".".inp" );

          $RLOC_file_to_run     = "$input"."_"."RGHT"."_"."$i";
	  $RLOC_out_file        = "$input"."_"."RGHT"."_"."$i".".log";
       system " rm ./$RLOC_file_to_run.dat"; 

	print GAMESSINPFILE2 "! \n";
	print GAMESSINPFILE2 "! \n";
	print GAMESSINPFILE2 "! \n";
	print GAMESSINPFILE2 "! \n";
	print GAMESSINPFILE2 "! \n";
	print GAMESSINPFILE2 "! \n";
	print GAMESSINPFILE2 "! \n";
        print GAMESSINPFILE2 " \$CONTRL SCFTYP=ROHF RUNTYP=ENERGY MAXIT=200         \n";
        print GAMESSINPFILE2 "          INTTYP=HONDO ICHARG=0 MULT=2 \$END\         \n";
        print GAMESSINPFILE2 " \$SYSTEM TIMLIM=3000 MEMORY=60000000  \$END\         \n";
        print GAMESSINPFILE2 " \$INTGRL SCHWRZ=.F. \$END\                           \n";
        print GAMESSINPFILE2 " \$NEO    NEOSCF=DIAGZN BASNUC=QZSPDD NUNIQN=1        \n";
	print GAMESSINPFILE2 "          IUNIQN(1)=25  EXCH=.T. NUMULT=2 LOCORB=2    \n";
        print GAMESSINPFILE2 "          NTAUXB=1 NAUXNB=1 IAUXNB(1)=26 \$END\       \n";
	print GAMESSINPFILE2 " \$BASIS  GBASIS=N31 NGAUSS=6 \$END\                  \n";
 	print GAMESSINPFILE2 " \$GUESS  GUESS=MOREAD NORB=150 \$END\                \n";
	print GAMESSINPFILE2 " \$DATA\ \n";
	print GAMESSINPFILE2 " \n";
        print GAMESSINPFILE2 " C1 \n";
        print GAMESSINPFILE2 " O 8.0    0.000000     0.000000    -1.200748 \n"; 
        print GAMESSINPFILE2 " O 8.0    0.000000     0.000000     1.200748 \n"; 
	print GAMESSINPFILE2 " C 6.0   -1.114666    -0.400558     1.749899 \n"; 
	print GAMESSINPFILE2 " C 6.0   -1.307180    -0.207036     3.150282 \n"; 
	print GAMESSINPFILE2 " C 6.0   -2.479880    -0.616856     3.759596 \n"; 
	print GAMESSINPFILE2 " C 6.0   -3.498257    -1.226627     3.007797 \n"; 
	print GAMESSINPFILE2 " C 6.0   -3.323431    -1.426760     1.629474 \n"; 
	print GAMESSINPFILE2 " C 6.0   -2.153726    -1.033298     1.000684 \n"; 
	print GAMESSINPFILE2 " C 6.0    1.114666    -0.400558    -1.749899 \n"; 
	print GAMESSINPFILE2 " C 6.0    2.153727    -1.033309    -1.000695 \n"; 
	print GAMESSINPFILE2 " C 6.0    3.323415    -1.426764    -1.629516 \n"; 
	print GAMESSINPFILE2 " C 6.0    3.498237    -1.226624    -3.007839 \n"; 
	print GAMESSINPFILE2 " C 6.0    2.479865    -0.616833    -3.759629 \n"; 
	print GAMESSINPFILE2 " C 6.0    1.307176    -0.207021    -3.150281 \n"; 
	print GAMESSINPFILE2 " H 1.0   -0.510170     0.275234     3.707418 \n"; 
	print GAMESSINPFILE2 " H 1.0   -2.616289    -0.462710     4.826741 \n"; 
	print GAMESSINPFILE2 " H 1.0   -4.415323    -1.547947     3.492806 \n"; 
	print GAMESSINPFILE2 " H 1.0   -4.109267    -1.904896     1.050621 \n"; 
	print GAMESSINPFILE2 " H 1.0   -2.002195    -1.195765    -0.061871 \n"; 
	print GAMESSINPFILE2 " H 1.0    2.002199    -1.195787     0.061859 \n"; 
	print GAMESSINPFILE2 " H 1.0    4.109246    -1.904915    -1.050671 \n"; 
	print GAMESSINPFILE2 " H 1.0    4.415300    -1.547945    -3.492852 \n"; 
	print GAMESSINPFILE2 " H 1.0    2.616271    -0.462682    -4.826774 \n"; 
	print GAMESSINPFILE2 " H 1.0    0.510165     0.275257    -3.707408 \n"; 
	print GAMESSINPFILE2 " H 1.0    0.000000     0.000000     $zh2     \n"; 
	print GAMESSINPFILE2 " H 1.0    0.000000     0.000000     $zh      \n"; 
	print GAMESSINPFILE2 " \$END\ \n";

        for( $k = 1; $k <= $RowsInOF2; ++$k ) {

        print GAMESSINPFILE2 "$RowsFromOrbitalFile2{$k}"; }

        $command_2 = " $GMS $RLOC_file_to_run $VER >& $RLOC_out_file";
        system( $command_2 );

#
#   Create the NEO-NCI input file
#

      open( GAMESSINPFILE3, ">$input"."_"."NOCI"."_"."$i".".inp" );

        $NOCI_file_to_run     = "$input"."_"."NOCI"."_"."$i";
        $NOCI_out_file        = "$input"."_"."NOCI"."_"."$i".".log";
        system " rm ./$NOCI_file_to_run.dat"; 

	print GAMESSINPFILE3 "! \n";
	print GAMESSINPFILE3 "! \n";
	print GAMESSINPFILE3 "! \n";
	print GAMESSINPFILE3 "! \n";
	print GAMESSINPFILE3 "! \n";
	print GAMESSINPFILE3 "! \n";
	print GAMESSINPFILE3 "! \n";
        print GAMESSINPFILE3 " \$CONTRL SCFTYP=ROHF RUNTYP=ENERGY MAXIT=200       \n";
        print GAMESSINPFILE3 "          INTTYP=HONDO ICHARG=0 MULT=2 \$END\       \n";
        print GAMESSINPFILE3 " \$SYSTEM TIMLIM=3000 MEMORY=125000000  \$END\      \n";
        print GAMESSINPFILE3 " \$INTGRL SCHWRZ=.F. \$END\                         \n";
        print GAMESSINPFILE3 " \$NEO    NEOSCF=DIAGZN BASNUC=QZSPDD NUNIQN=1      \n";
	print GAMESSINPFILE3 "          NTAUXB=1 NAUXNB=1 NUMULT=2 NEONCI=.T.     \n";
	print GAMESSINPFILE3 "          IUNIQN(1)=25 IAUXNB(1)=26 EXCH=.T. \$END\ \n";
	print GAMESSINPFILE3 " \$NCIINP NCIPRT=3   \$END\                         \n";
	print GAMESSINPFILE3 " \$BASIS  GBASIS=N31 NGAUSS=6 \$END\                \n";
 	print GAMESSINPFILE3 " \$GUESS  GUESS=MOREAD NORB=150 \$END\              \n";
	print GAMESSINPFILE3 " \$DATA\ \n";
	print GAMESSINPFILE3 " \n";
        print GAMESSINPFILE3 " C1                                              \n";
        print GAMESSINPFILE3 " O 8.0    0.000000     0.000000    -1.200748 \n"; 
        print GAMESSINPFILE3 " O 8.0    0.000000     0.000000     1.200748 \n"; 
	print GAMESSINPFILE3 " C 6.0   -1.114666    -0.400558     1.749899 \n"; 
	print GAMESSINPFILE3 " C 6.0   -1.307180    -0.207036     3.150282 \n"; 
	print GAMESSINPFILE3 " C 6.0   -2.479880    -0.616856     3.759596 \n"; 
	print GAMESSINPFILE3 " C 6.0   -3.498257    -1.226627     3.007797 \n"; 
	print GAMESSINPFILE3 " C 6.0   -3.323431    -1.426760     1.629474 \n"; 
	print GAMESSINPFILE3 " C 6.0   -2.153726    -1.033298     1.000684 \n"; 
	print GAMESSINPFILE3 " C 6.0    1.114666    -0.400558    -1.749899 \n"; 
	print GAMESSINPFILE3 " C 6.0    2.153727    -1.033309    -1.000695 \n"; 
	print GAMESSINPFILE3 " C 6.0    3.323415    -1.426764    -1.629516 \n"; 
	print GAMESSINPFILE3 " C 6.0    3.498237    -1.226624    -3.007839 \n"; 
	print GAMESSINPFILE3 " C 6.0    2.479865    -0.616833    -3.759629 \n"; 
	print GAMESSINPFILE3 " C 6.0    1.307176    -0.207021    -3.150281 \n"; 
	print GAMESSINPFILE3 " H 1.0   -0.510170     0.275234     3.707418 \n"; 
	print GAMESSINPFILE3 " H 1.0   -2.616289    -0.462710     4.826741 \n"; 
	print GAMESSINPFILE3 " H 1.0   -4.415323    -1.547947     3.492806 \n"; 
	print GAMESSINPFILE3 " H 1.0   -4.109267    -1.904896     1.050621 \n"; 
	print GAMESSINPFILE3 " H 1.0   -2.002195    -1.195765    -0.061871 \n"; 
	print GAMESSINPFILE3 " H 1.0    2.002199    -1.195787     0.061859 \n"; 
	print GAMESSINPFILE3 " H 1.0    4.109246    -1.904915    -1.050671 \n"; 
	print GAMESSINPFILE3 " H 1.0    4.415300    -1.547945    -3.492852 \n"; 
	print GAMESSINPFILE3 " H 1.0    2.616271    -0.462682    -4.826774 \n"; 
	print GAMESSINPFILE3 " H 1.0    0.510165     0.275257    -3.707408 \n"; 
	print GAMESSINPFILE3 " H 1.0    0.000000     0.000000     $zh2     \n"; 
	print GAMESSINPFILE3 " H 1.0    0.000000     0.000000     $zh      \n"; 
	print GAMESSINPFILE3 " \$END\ \n";

        for( $k = 1; $k <= $RowsInOF; ++$k ) {
                                                                          
        print GAMESSINPFILE3 "$RowsFromOrbitalFile{$k}"; }
# 
#     Get the left localized orbitals from the DAT file
#
         $stringlefte  = "LEFT LOCALIZED ELEC"; 
         $stringleftp  = "LEFT LOCALIZED NUC"; 
         $stringrighte = "RIGHT LOCALIZED ELEC"; 
         $stringrightp = "RIGHT LOCALIZED NUC"; 
         $stringend    = "\$END"; 
         $switch1 = 0;
         $stop1   = 0;         
         open( FILETOREAD, "<$input"."_"."LEFT"."_"."$i".".dat" );
	 while( $line = <FILETOREAD> ) {
		 
		if(  $line =~ /$stringlefte/ || $line =~ /$stringleftp/ ) {
                       $switch1 = 1;         } 
                
		if( $switch1 == 1 )            {
                   print GAMESSINPFILE3 $line;
		if( $line =~ /$stringend/ ) {
                       $switch1 = 0;         } }

                  ++$line;             }

# 
#     Get the right localized orbitals from the DAT file
#
         $switch1 = 0;
         open( FILETOREAD, "<$input"."_"."RGHT"."_"."$i".".dat" );
	 while( $line = <FILETOREAD> ) {
		 
		if(  $line =~ /$stringrighte/ || $line =~ /$stringrightp/ ) {
                       $switch1 = 1;         } 
                
		if( $switch1 == 1 )            {
                   print GAMESSINPFILE3 $line;
		if( $line =~ /$stringend/ ) {
                       $switch1 = 0;         } }

                  ++$line;             }
#
#     Execute the NEO-NCI job
#

         $command_3 = " $GMS $NOCI_file_to_run $VER >& $NOCI_out_file";
      
         system( $command_3 );
#
#     Assign current H basis function center separation to array dist[]
#

   $dist[$i] = -$zh2 + $zh;

# 
# GET ALL RELEVANT ENERGIES AND SPLITTINGS
#

         open( FILETOREAD, "<$input"."_"."NOCI"."_"."$i".".log" );

	  while( $line = <FILETOREAD> ) {
		 
		if( $line =~ /$string1/ ) {
			@line_chunks = split( /\s+/, $line );
                    $es1[$i] = $line_chunks[3]; 
	            $es2[$i] = $line_chunks[4];
	            $h11[$i] = $line_chunks[5];
	            $h12[$i] = $line_chunks[6]; 
	            $S12[$i] = $line_chunks[7]; }
                ++$line; } 

#
# GET THE VIBRONIC SPLITTING
#
	     open(FILETOREAD, "<$input"."_NOCI_"."$i".".log" );

	     while( $line = <FILETOREAD> ) {
               if( $line =~ /$string2/ ) {
                    @line_chunks2 = split( /\s+/, $line );
                    $DE[$i] = $line_chunks2[4] ;}
                   ++$line; } 

   printf MATH_DATAFILE "    %6.4f\t  %14.8f\t  %14.8f\t %14.8f\t %14.8f\t %10.6f\t %9.3f \
                         \n", $dist[$i], $es1[$i], $es2[$i], $h11[$i], $h12[$i], $S12[$i], $DE[$i]; 

#
#  GET THE NONORTHOGONAL OVERLAPS


        open( FILETOREAD, "<$input"."_"."NOCI"."_"."$i".".log" );

	 while( $line = <FILETOREAD> ) {
		 
             if( $line =~ /$string3/ ) {
	         @line_chunks3 = split( /\s+/, $line );
                $Se12[$i] = $line_chunks3[4]; 
	         $Sp12[$i] = $line_chunks3[5];}

                ++$line; } 

   print MATH_DATAFILE2 "       $dist[$i]\t   $Se12[$i]\t   $Sp12[$i]\t  $S12[$i]\t  $DE[$i]\n"; 

#
#  Update the proton coordinates
#
      $zh  = $zh  + $disph1;
      $zh2 = $zh2 + $disph2; 
      }
 }

 exit
