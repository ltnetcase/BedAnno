use 5.010001;
use ExtUtils::MakeMaker;

# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
	NAME          => 'BedAnno',
	VERSION_FROM  => 'lib/BedAnno.pm',            # finds $VERSION
	PREREQ_PM     => { 'Bio::DB::HTS' => '>=3.00', 'Time::HiRes' => '>=1.90', 'IO::Uncompress::Gunzip' => '>=2.00'},    # e.g., Module::Name => 1.1
	TEST_REQUIRES => { 'Test::Most' => '>=0.34' },    # for test only
	(
		$] >= 5.005
		?    ## Add these new keywords supported since 5.005
		  (
			ABSTRACT_FROM => 'lib/BedAnno.pm',   # retrieve abstract from module
			AUTHOR => 'liutao <ltnetcase@gmail.com>'
		  )
		: ()
	),
);
