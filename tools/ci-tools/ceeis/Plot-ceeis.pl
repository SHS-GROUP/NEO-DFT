#!/usr/bin/perl
use File::Basename;
use Getopt::Long;


$DoPlot=1;
$DoMatrix=1;
$DoTriples=0;
$DoPostscript=0;

usage() if ( ! GetOptions('help|?' => \$help,
		'plot!' => \$DoPlot,
		'gnuplot-file=s' => \$GnuplotFile,
		'matrix!' => \$DoMatrix,
		'matrix-file=s' => \$MatrixFile,
		'triples!' => \$DoTriples,
		'postscript!' => \$DoPostscript)

		or @ARGV < 1 or defined $help);


#this will be used by &PrintMatrix and &PrintGnuplot
#to determine file names.
($GAMESSout, $path, $suffix) = fileparse( @ARGV[0], qr/\.[^.]*/ );


#Some global parameters that will 
#be usefull for line searching.
$HighestExcitation;
$HighestOrbital;
$LowestVirtual;

#this is necessary when using
#GetExcitations
$MatrixCount=1;


#this will contain DEMAT(m,x) this is a hashtable beacuase
#larger excitations may not contain as many elements as 
#the double and tripple and I don't want to plot 0.0
%Excitations=();

#For each relevant excitation, the fitting parameters a*x+b 
#will be stored.  If fitting paramters are not found, gnuplot
#will optimize them.  The fitting iterations will be printed
#to '$GnuplotFile.log'.  There may be a way to have gnuplot write this
#to another name.
%Fitting=(); 

open(in, @ARGV[0]);
while(<in>)
{


	#Once we have the Global variables, 
	#we can go find stuff.
	&GetGlobals($_);

	#Set the keys in %Excitations.
	&PrepareExcitations($_);

	#Set the values in %Excitations.
	&GetExcitations($_);

	#Get extrapolation parameters, if they exist.
	&GetFitting($_);


}
close(in);

unless ($DoMatrix==0)
{
	&PrintMatrix;
}

unless ($DoPlot==0)
{
	&PrintGnuplot;
}




#------------------------------ SUBROUTINES -----------------------------------------


sub GetGlobals
{
	#The following if statements set the necessary global variables

	#NORB is the size of the basis set and incidently the orbital number
	#to which we are extrapolating.
	if($_[0]=~/NORB  =/)
	{
		@line=split("NORB  =",$_[0]);
		$temp=@line[1];
		@line2=split(" ",$temp);
		$HighestOrbital=@line2[0];
		print STDOUT "Highest oribtal is $HighestOrbital\n";
	}


	#MMIN is the lowest virtual orbital that will appear in the DEMAT
	if($_[0]=~/MMIN  =/)
	{
		@line=split("MMIN  =",$_[0]);
		$temp=@line[1];
		@line2=split(" ",$temp);
		$LowestVirtual=@line2[0];
		print STDOUT "Lowest virtual orbital is $LowestVirtual\n";
	}

	#ISTPEX is the highest excitation which will be extrapolated.
	#Practically, only excitations 3,4,5,6 can be easily plotted on
	#one page.
	if($_[0]=~/ISTPEX=/)
	{
		@line=split("ISTPEX=",$_[0]);
		$temp=@line[1];
		@line2=split(" ",$temp);
		$HighestExcitation=@line2[0];
		print STDOUT "Highest excitation is $HighestExcitation\n";
	}
}

sub PrepareExcitations
{

	#The explicitly computed elements of DEMAT will be placed in the
	#hashtable %Excitations{$ex}{$orb}, where $ex is 2, 3, 4, etc.
	#and $orb is the correlating orbital number.  For each excitation,
	#the $HighestVirtual is also included even though this is extrapolated
	#for $ex >= 3.

	#prepare the Excitations which will keep the running total
	#of excitation energies for each level as specified in M1M2ex
	if($_[0]=~ /EXCITATION LEVELS     HIGHEST ORBITALS USED/)
	#After this line, the orbitals used for extrapolation can appear
	#as a discrete list or upper and lower bounds.  Either way, we 
	#are going to need the discrete list.  While, each excitation can
	#be computed with any number of orbitals, it is only those that match
	#the X-2 excitation that will be plotted.
	{
		#Read ahead for X= 
		%Hold=();
		$Hold{2}=1; #ex=2 has been computed for all orbitals
		$Highest=0; #This ensures we don't run to the end of file
		do {
			$nextline=<in>;
			if($nextline=~/X=/)
			{
				@line=split(" ",$nextline);
				$ex=@line[1];
				#Have we reached the ISTPEX excitation?
				if($ex == $HighestExcitation) {
					$Highest=1;
				}
				@line2=split("M=",$nextline);
				$orbs=@line2[1];
				@orbs=split(" ",$orbs);
				#Previous is used when orbitals appear
				#as bounds.
				$previous=-1;
				chop($orbs);
				print STDOUT "I found excitation level: $ex\n";
				print STDOUT " with orbs $orbs\n";
				#we musn't forget that $HighestOrbital is important too
				#TODO: when orbs= 0,0
				foreach $i (@orbs,$HighestOrbital)
				{
					if($i==0 && $previous==0)
					{
						#this is the case that all orbs are computed
						#we will populate these as we go
						print STDOUT "I will hold level: $ex\n";
						$Hold{$ex}=1;
					}
					elsif($i==0 && $previous!=0)
					{
						$previous=$i;
					}
					elsif($i<0)
					#Exapand the orbital numbers.
					{
						foreach $j ($previous+1 .. -$i)
						{
							$Excitations{$ex}{$j}=0.0;
							foreach $hold (keys %Hold)
							{
								$Excitations{$hold}{$j}=0.0;
							}
						}

					}
					elsif($i>0)
					{
						$Excitations{$ex}{$i}=0.0;
						$previous=$i;
						foreach $hold (keys %Hold)
						{
							$Excitations{$hold}{$i}=0.0;
						}
					}

				}
			}
		} while($Highest==0);
	} 

}

sub GetExcitations
{

	#$HighestExcitation-1 Matrices will appear in the output.
	#The first one is excitation level 2 and will only contain
	#the items necessary to compare to higher excitations.
	if($_[0]=~/CEEIS DIFFERENCE MATRIX DEMAT/)
	{

		#Because DEMAT will appear $HighestExcitation-1 times without
		#announcing nearby which excitation level was most recently
		#calculated, $MatrixCount will be used to set values in
		#%Excitations.
		$MatrixCount++;


		print "Populating computed excitations in level $MatrixCount\n";
		#Starting from $LowestVirtual, populate excitation energies
		#only for the orbitals in M1M2ex.
		$nextline=<in>;


		#In older versions of GAMESS,
		#the format statement was set 
		#to I2, so a second counter is
		#necessary when the number of 
		#orbitals exceeds 99.
		#Since late 2008, the format statement
		#has been changed to I4.  
		$count=$LowestVirtual;
		foreach $i($LowestVirtual .. $HighestOrbital)
		{
			$nextline=<in>;
			@line=split(" ",$nextline);
			$orbital=@line[0];
			$energy=@line[$MatrixCount-1];

			foreach $i (sort keys %{$Excitations{$MatrixCount}})
			{
				if($orbital=~/\*\*/ && $count==$i)
				{
					$Excitations{$MatrixCount}{$i}=$energy;
				}

				if($i==$orbital)
				{
					$Excitations{$MatrixCount}{$i}=$energy;
				}
			}


			$count++; 

		}

	}

}

sub GetFitting
{
	if($_[0]=~/EXTRAPOLATION FOR EXCITATION LEVEL/)
	{
		@line=split(" ",$_[0]);
		$ex=@line[4];

		$finished=0;
		do
		{
			$nextline=<in>;
			if($nextline=~/Y =/)
			{
				@line=split(" ",$nextline);
				$a=@line[2];
				@line2=split("X",$a);
				$a=@line2[0];
				$b=@line[4];

				$Fitting{$ex}{'a'}=$a;
				$Fitting{$ex}{'b'}=$b;
				$finished=1;

				print "Excitation $ex has been fit with Y = $a*x+$b\n";
			}
		} while($finished==0);
	}
}

sub PrintMatrix
{


	#How shall we name the file?
	if(!$MatrixFile)
	{
		$MatrixFile=$GAMESSout.".txt";
		print STDOUT "Matrix file will be written to $MatrixFile\n";
	}
	else
	{
		print STDOUT "Matrix file will be written to $MatrixFile\n";
	}

	#Print the Excitations to the data file for 
	#plottting to $MatrixFile.
	#Ex=2 must always exist and contain all possible 
	#orbitals for use in other excitations.  If an
	#orbital does not exist, nothing will be printed.

	open(matrix,">$MatrixFile");
	print matrix "#Orbital      ";
	foreach $ex (sort {$a<=>$b} keys %Excitations)
	{
		print matrix "Ex: $ex            ";
	}
	print matrix "\n";
	foreach $orb (sort {$a<=>$b} keys %{$Excitations{2}})
	{
		printf matrix "%-4d     ",$orb;

		foreach $ex (sort {$a<=>$b} keys %Excitations)
		{
			printf matrix "%10.12f  ",$Excitations{$ex}{$orb};

		}

		print matrix "\n";
	}
	close(matrix);

}


sub PrintGnuplot
{

	#Oh, dear.  These names can get complicated.
	if(!$GnuplotFile)
	{
		$GnuplotFile=$GAMESSout.".gnuplot";
		print STDOUT "GNUPLOT file will be written to $GnuplotFile\n";
	}
	else
	{
		print STDOUT "GNUPLOT file will be written to $GnuplotFile\n";
	}
	if($DoPostscript)
	{
		$PostscriptFile=$GAMESSout.".ps";
	}



	#print the gnuplot file
	open(gnuplot,">$GnuplotFile");
	if($DoPostscript)
	{
		print gnuplot
		"set term postscript enhanced solid color
set out '$PostscriptFile'\n"
	}
	print gnuplot 
	"set size 1,1
set zero 1e-20
set origin 0,0
set multiplot\n";

	#The following is to plot Triple Excitations against Double Excitations.
	if($DoTriples==1)
	{
		print gnuplot "set size 0.5,0.5\nset origin 0,0.5\n";
		print gnuplot "set title 'E(3|m) : E(2|m)'\n";
		print gnuplot "set xlabel '-{/Symbol D}E(2|m) mhartree'\n";
		print gnuplot "set ylabel '-{/Symbol D}E(3|m) mhartree'\n";
		if( exists $Fitting{3})
		{
			printf gnuplot "plot '%s' u (-\$2*1000.0):(-\$3*1000.0) notitle, %e*x-%e*1000.0 notitle\n",$MatrixFile,$Fitting{3}{a},$Fitting{3}{b};
		}
		else
		{
			#Triples are not always fitted against Doubles,
			#so GNUPlot must do it.
			print gnuplot "a3=1.0\n";
			print gnuplot "b3=1.0\n";
			print gnuplot "f3(x)=a3*x+b3\n";
			print gnuplot "fit f3(x) '$MatrixFile' using 2:3 via a3,b3\n";
			print gnuplot "plot '$MatrixFile' u (-\$2*1000.0):(-\$3*1000.0) notitle, a3*x-1000.0*b3 notitle\n";
		}
	}

	if (exists $Fitting{4})
	{
		print gnuplot "set size 0.5,0.5\n";
		if ($DoTriples==1) 
		#How many plots are we going to make?
		{
			print gnuplot "set origin 0.5,0.5\n";
		}
		else
		{
			print gnuplot "set origin 0.25,0.5\n";
		}
		print gnuplot "set xlabel '-{/Symbol D}E(2|m) mhartree'\n";
		print gnuplot "set ylabel '-{/Symbol D}E(4|m) mhartree'\n";
		print gnuplot "set title 'E(4|m) : E(2|m)'\n";
		printf gnuplot "plot '%s' u (-\$2*1000.0):(-\$4*1000.0) notitle, %e*x-%e*1000.0 notitle\n",$MatrixFile,$Fitting{4}{a},$Fitting{4}{b};
	}

	if(exists $Fitting{5})
	{
		print gnuplot "set size 0.5,0.5\nset origin 0,0\n";
		print gnuplot "set title 'E(5|m) : E(3|m)'\n";
		print gnuplot "set xlabel '-{/Symbol D}E(3|m) mhartree'\n";
		print gnuplot "set ylabel '-{/Symbol D}E(5|m) mhartree'\n";
		printf gnuplot "plot '%s' u (-\$3*1000.0):(-\$5*1000.0) notitle, %e*x-%e*1000.0 notitle\n",$MatrixFile,$Fitting{5}{a},$Fitting{5}{b};
	}

	if(exists $Fitting{6})
	{
		print gnuplot "set size 0.5,0.5\nset origin 0.5,0\n";
		print gnuplot "set title 'E(6|m) : E(4|m)'\n";
		print gnuplot "set xlabel '-{/Symbol D}E(4|m) mhartree'\n";
		print gnuplot "set ylabel '-{/Symbol D}E(6|m) mhartree'\n";
		printf gnuplot "plot '%s' u (-\$4*1000.0):(-\$6*1000.0) notitle, %e*x-%e*1000.0 notitle\n",$MatrixFile,$Fitting{6}{a},$Fitting{6}{b};
	}

	print gnuplot "unset multiplot\n";

	close(gnuplot);

}


sub usage
{
	print STDOUT "
	usage: plot-ceeis [options] gamess-output

	Options
	--no-plot                   Disable writting of GNUPlot file.
	--gnuplot-file=filename     Name and location of GNUPlot file.
	                            The default is gamess-output.gnuplot.
	--no-matrix                 Disable writting of Difference Matrix file.
	--matrix-file=filename      Name and location of Difference Matrix file.
	                            The default is gamess-output.txt.
	--triples                   Plot Triple Excitations against Double Excitations.
	--postscript                Print GNUPlot directly to postcript file.


	Albert DeFusco 2008 Ames Laboratory\n";
	exit;
}
