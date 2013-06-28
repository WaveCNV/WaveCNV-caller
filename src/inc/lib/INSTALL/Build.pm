#------------------------------------------------------------------------
#----                         INSTALL::Build                         ----
#------------------------------------------------------------------------
package INSTALL::Build;
use strict;
use warnings;
use POSIX;
use Config;
use DynaLoader;
use File::Copy;
use File::Path;
use File::Which;
use vars qw($BIN);
use Cwd ();

BEGIN{
    #prepare correct version of Module Build for building everything
    my $Bundled_MB = 0.3607;  #version included in my distribution

    # Find out what version of Module::Build is installed or fail quietly.
    # This should be cross-platform.
    my $Installed_MB =`$^X -e "eval q{require Module::Build; print Module::Build->VERSION} or exit 1"`;
    chomp $Installed_MB;
    $Installed_MB = 0 if $?;

    $BIN = Cwd::cwd();

    # Use the bundled copy of Module::Build if it's newer than the installed.
    if ($Bundled_MB > $Installed_MB){
	unshift @INC, "$BIN/inc/bundle" unless($INC[0] eq "$BIN/inc/bundle");
	my $PERL5LIB = $ENV{PERL5LIB} || '';
	$PERL5LIB = "$BIN/inc/bundle:$PERL5LIB";
	$ENV{PERL5LIB} = $PERL5LIB;
    }

    require Module::Build;
}

use base qw(Module::Build);
__PACKAGE__->add_property( 'exe_requires' );
__PACKAGE__->add_property( 'lib_requires' );

eval 'require LWP::Simple';
eval 'require Archive::Tar';

#------------------------------------------------------------------------
#--------------------------------- METHODS ------------------------------
#------------------------------------------------------------------------
sub new {
    my $class = shift @_;

    my %hash = @_;
    my $self = $class->SUPER::new(@_);

    #fix weird override bug on systems with local::lib
    if($hash{install_base}){
	$self->install_base($hash{install_base});
	$self->config_data(install_base => $hash{install_base});
    }
    elsif($self->config_data('install_base')){
	$self->install_base($self->config_data('install_base'));
    }

    $self->install_base_relpaths('exe'  => 'exe');
    $self->install_base_relpaths('data' => 'data');
    bless($self, $class);

    #performs a check for eternal algorithm dependencies
    $self->check_exes;
    $self->check_libs;

    return $self;
}

sub resume {
    my $class = shift @_;
    my $self = $class->SUPER::resume(@_);

    #fix weird override bug on systems with local::lib
    if($self->config_data('install_base')){
	$self->install_base($self->config_data('install_base'));
    }

    return $self;
}

#returns MPI compiler and includes directory
sub config_mpi {
    my $self = shift;

    my $base = $self->base_dir;
    my $ebase = $self->install_destination('exe');
    my @exes = grep {/(^|[\/])mpicc$/} (<$base/../exe/*/*>, <$base/../exe/*/bin/*>);
    my $mpicc = "$ebase/mpich2/bin/mpicc" if(-f "$ebase/mpich2/bin/mpicc");
    ($mpicc) = grep {/(^|[\/])mpicc$/} (<$base/../exe/*/*>, <$base/../exe/*/bin/*>) if(!$mpicc);
    $mpicc = $self->config('cc') if(! $mpicc && $self->config('cc') =~ /(^|[\/])mpicc$/);
    ($mpicc) = File::Which::where('mpicc') if(!$mpicc || ! -f $mpicc);

    $mpicc = $self->prompt("\nPlease specify the path to 'mpicc' on your system:", $mpicc);

    while($mpicc !~ /(^|[\/])mpicc$/ || ! -f $mpicc){
	$mpicc = $self->prompt("\nCannot find 'mpicc'.\n".
			    "Please specify the path (leave blank to skip):", '');
	return if(! $mpicc);
    }

    my $ccdir = $mpicc;
    $ccdir =~ s/\/+[^\/]+\/[^\/]+$//;

    #directories to search for mpi.h
    my @includes = (<$ccdir/include>,
		    <$ebase/mpich2/include>,
		    <$ebase/*/include>,
		    </usr/include>,
		    </usr/include/mpi*>,
		    </usr/mpi*/include>,
		    </usr/local/include>,
		    </usr/local/include/mpi*>,
		    </usr/local/mpi*/include>,
		    </usr/lib/>,
		    </usr/lib/include/mpi*>,
		    </usr/lib/mpi*/include>,
		    </usr/local/lib>,
		    </usr/local/lib/include/mpi*>,
		    </usr/local/lib/mpi*/include>);

    my ($MPIDIR) = grep {-f "$_/mpi.h"} @includes;

    $MPIDIR = $self->prompt("\nPlease specify the path to the directory containing 'mpi.h':", $MPIDIR);

    while(! -f "$MPIDIR/mpi.h"){
	$MPIDIR = $self->prompt("\nCannot find 'mpi.h'\n.".
				"Please specify the containing directory path (leave blank to cancel):", '');
	return if(! $MPIDIR);
    }

    $self->config_data(MPIDIR => $MPIDIR);
    $self->config_data(MPICC => $mpicc);

    $self->add_exe_requires(mpicc => $mpicc);
    $self->add_lib_requires(MPI => "$MPIDIR/mpi.h");
}

#sets locations for SAM compiler files
sub config_sam_libs {
    my $self = shift;

    my $samdir = $self->config_data('samtools') || '';
    $samdir =~ s/\/[^\/]+$//;

    return if(!$samdir);

    #directories to search for bam.h and libbam.a
    $ENV{SAMTOOLS} ||= $samdir;
    my $ebase = $self->install_destination('exe');
    my @includes = (<$samdir/>,
		    <$samdir/include>,
		    <$samdir/lib>,
		    <$ebase/samtools/include>,
		    <$ebase/samtools/lib>,
		    <$ebase/*/include>,
		    <$ebase/*/lib>,
		    <$ENV{SAMTOOLS}/>,
		    </usr/include>,
		    </usr/include/sam*>,
		    </usr/include/lib>,
		    </usr/include/lib/sam*>,
		    </usr/sam*>,
		    </usr/sam*/include>,
		    </usr/sam*/lib>,
		    </usr/local/sam*>,
		    </usr/local/include>,
		    </usr/local/include/sam*>,
		    </usr/local/include/sam*/lib>,
		    </usr/local/sam*/include>,
		    </usr/local/lib>,
		    </usr/local/lib/sam*>,
		    </usr/local/sam*/lib>,
		    </usr/share>,
		    </usr/share/sam*>,
		    </usr/share/include>,
		    </usr/share/include/sam*>,
		    </usr/share/sam*/include>,
		    </usr/share/lib>,
		    </usr/share/lib/sam*>,
		    </usr/share/sam*/lib>,
		    </usr/lib/>,
		    </usr/lib/sam*>,
		    </usr/lib/include>,
		    </usr/lib/include/sam*>,
		    </usr/lib/sam*/include>,
		    </usr/local/lib>,
		    </usr/local/lib/sam*>,
		    </usr/local/lib/include/sam*>,
		    </usr/local/lib/sam*/include>);

    my ($sam_include) = grep {-f "$_/bam.h"} @includes;
    $sam_include = ($sam_include) ? "$sam_include/bam.h" : '';
    $sam_include =~ s/\/+/\//g;
    $sam_include = $self->prompt("\nPlease specify the path to samtools 'bam.h' on your system:", $sam_include);
    return if(!$sam_include);
    while(! -f $sam_include){
	$sam_include = $self->prompt("\nCannot find 'bam.h'\n.".
				"Please specify the path (leave blank to skip):", '');
	return if(!$sam_include);
    }

    my ($sam_lib) = grep {-f "$_/libbam.a"} @includes;
    $sam_lib = ($sam_lib) ? "$sam_lib/libbam.a" : '';
    $sam_lib =~ s/\/+/\//g;
    $sam_lib = $self->prompt("\nPlease specify the path to samtools 'libbam.a' on your system:", $sam_lib);
    return if(!$sam_lib);
    while(! -f $sam_lib){
	$sam_lib = $self->prompt("\nCannot find 'libbam.a'\n.".
				"Please specify the path (leave blank to skip):", '');
	return if(! $sam_lib);
    }

    $sam_include =~ s/[^\/]+$//; #make path to containing directory
    $sam_lib =~ s/[^\/]+$//; #make path to containing directory
    $self->config_data(sam_include => $sam_include);
    $self->config_data(sam_lib => $sam_lib);
    $self->config_data(include_extra => [$sam_include]);
}

sub config_exe_loc {
    my $self = shift;
    my $tag = shift;

    my $base = $self->base_dir;
    my $ebase = $self->install_destination('exe');
    my ($exe) = grep {/(^|[\/])$tag$/ && -f $_} (<$base/../exe/*/*>, <$base/../exe/*/bin/*>);
    ($exe) = File::Which::where($tag) if(!$exe || ! -f $exe);
    $exe = $self->prompt("\nPlease specify the path to '$tag' on your system:", $exe);
    if(!$exe){
	print "Skipping '$tag'...\n";
	return $exe;
    }
    while($exe !~ /(^|[\/])$tag$/ || ! -f $exe){
        $exe = $self->prompt("\nCannot find '$tag'.\n".
			     "Please specify the path (leave blank to skip):", '');
        return if(! $exe);
    }

    $self->config_data($tag => $exe);
    return $exe;
}



#add a Module to the requires list
sub add_requires {
    my $self = shift;
    return unless(@_);

    if (@_ > 1 && @_ % 2 == 0){
	for (my $i = 0; $i < @_; $i += 2){
	    $self->{properties}{requires}{$_[$i]} = $_[$i+1]
	}
    }
    elsif(ref $_[0] eq 'HASH'){
	map {$self->{properties}{requires}{$_} = $_[0]->{$_} } keys %{$_[0]}
    }
    else{
	map {$self->{properties}{requires}{$_} = 0} @_;
    }
}

#add a Module to the build_requires list
sub add_build_requires {
    my $self = shift;
    return unless(@_);

    if (@_ > 1 && @_ % 2 == 0){
	for (my $i = 0; $i < @_; $i += 2){
	    $self->{properties}{build_requires}{$_[$i]} = $_[$i+1]
	}
    }
    elsif(ref $_[0] eq 'HASH'){
	map {$self->{properties}{build_requires}{$_} = $_[0]->{$_} } keys %{$_[0]}
    }
    else{
	map {$self->{properties}{build_requires}{$_} = 0} @_;
    }
}

#add a Module to the exe_requires list
sub add_exe_requires {
    my $self = shift;
    return unless(@_);

    if (@_ > 1 && @_ % 2 == 0){
	for (my $i = 0; $i < @_; $i += 2){
	    $self->{properties}{exe_requires}{$_[$i]} = $_[$i+1]
	}
    }
    elsif(ref $_[0] eq 'HASH'){
	map {$self->{properties}{exe_requires}{$_} = $_[0]->{$_} } keys %{$_[0]}
    }
    else{
	map {$self->{properties}{exe_requires}{$_} = 0} @_;
    }
}

#add a Library to the lib_requires list
sub add_lib_requires {
    my $self = shift;
    return unless(@_);

    if (@_ > 1 && @_ % 2 == 0){
	for (my $i = 0; $i < @_; $i += 2){
	    $self->{properties}{lib_requires}{$_[$i]} = $_[$i+1]
	}
    }
    elsif(ref $_[0] eq 'HASH'){
	map {$self->{properties}{lib_requires}{$_} = $_[0]->{$_} } keys %{$_[0]}
    }
    else{
	map {$self->{properties}{lib_requires}{$_} = 0} @_;
    }
}

#replaces Module::Build's config method
sub config {
    my $self = shift;
    
    #hack to get around bug in Module::Build 
    #otherwise it will ignore changes to config 'cc' adn 'ld'
    $self->{stash}->{_cbuilder} = undef;

    return $self->SUPER::config(@_);
}

#override default build
sub ACTION_build {
    my $self = shift;

    my @perl = map {keys %{$_->{requires}}} $self->prereq_failures();
    my @exes = map {keys %{$_->{exe_requires}}} $self->exe_failures();
    my @libs = map {keys %{$_->{lib_requires}}} $self->lib_failures();

    if(@perl || @exes || @libs){
	$self->status;
	print "\nERROR: Cannot '".$self->invoked_action."', missing prerequisites\n";
	exit;
    }
    
    if($self->feature('mpi_support')){
	$self->log_info("Configuring " . $self->dist_name . " with MPI support\n");
    }
    else{
	$self->log_info("Building " . $self->dist_name . "\n");
    }
    $self->depends_on('code');
    $self->depends_on('docs');

    #compile MPI module
    if($self->feature('mpi_support')){
        require Parallel::Application::MPI;
	my $loc = $self->blib();
	mkdir($loc) if(!$loc);
        Parallel::Application::MPI::_bind($self->config_data('MPICC'),
                                          $self->config_data('MPIDIR'),
					  $loc);	
	File::Path::rmtree($self->blib()."/build");
    }    

    $self->config_data(build_done => 1);
}

#commit current SRC to subversion repository
sub ACTION_commit {
    my $self = shift;

    $self->sync_bins();
    my ($s_svn) = `svn info` =~ /Revision\:\s*(\d+)/;
    $self->svn_w_args('update', '');
    my ($f_svn) = `svn info` =~ /Revision\:\s*(\d+)/;
    $self->svn_w_args('commit');

    #there were changes so re-run install
    if($s_svn != $f_svn){
	$self->dispatch('clean');
	$self->dispatch('install');
    }
}

#update current SRC from the SVN repository
sub ACTION_update {
    my $self = shift;

    $self->sync_bins();
    my ($s_svn) = `svn info` =~ /Revision\:\s*(\d+)/;
    $self->svn_w_args('update');
    my ($f_svn) = `svn info` =~ /Revision\:\s*(\d+)/;

    #there were changes so re-run install
    if($s_svn != $f_svn){
	$self->dispatch('clean');
	$self->dispatch('install');
    }

    print "\nSVN STATUS:\n";
    $self->svn_w_args('status');
}

#syncronize the .../src/bin and .../src/inc/bin directories
#to .../bin because of user edits to scripts
sub ACTION_sync {
    my $self = shift;

    $self->sync_bins();
}

#creates a versioned release
sub ACTION_release {
    my $self = shift;

    #update current SRC
    print "\nUpdating to most current repository...\n";
    $self->sync_bins();
    my ($s_svn) = `svn info` =~ /Revision\:\s*(\d+)/;
    $self->svn_w_args('update', '');

    #doing
    print "\nPre-release commit of any user changes...\n";
    $self->svn_w_args('commit', '-m "pre-release commit"');
    $self->svn_w_args('update', '');

    #clean and check versions for release
    $self->dispatch('clean');
    my $ver = $self->check_update_version(); #returns package version
    $self->{properties}->{dist_version} = $ver;

    File::Which::which('tar') || die "ERROR: Cannot find tar to build the release\n";
    File::Which::which('svn') || die "ERROR: Cannot find the executable svn\n";
    
    #build tarball for users to download
    my $cwd = $self->base_dir;
    my $tgz = "$cwd/package-$ver.tgz";    
    if(! -f $tgz){
	my ($dir, $base) = $cwd =~ /^(.*\/)([^\/]+)\/src$/;
	
	my $exclude = `cd $dir; svn status $base`;
	$exclude = join("\n", ($exclude =~ /\?\s+([^\n]+)/g), "src/package-$ver.tgz") ."\n";
	open(OUT, "> .exclude~");
	print OUT $exclude;
	close(OUT);
	
	print "\nBuilding tarball for distribution...\n";
	my $command = "tar -C $dir -zcf $tgz $base --exclude \"*~\" --exclude \".svn\" --exclude \"package-*.tgz\" --exclude-from .exclude~";
	system($command) && unlink($tgz);
	unlink(".exclude~");
	die "ERROR: tarball creation failed\n" if(! -f $tgz);
    }

    #there were changes so re-run install (updates version info in scripts)
    my ($f_svn) = `svn info` =~ /Revision\:\s*(\d+)/;
    if($s_svn != $f_svn){
	print "\nNow reinstalling scripts to reflect version changes...\n";
	sleep 1;
	$self->dispatch('realclean'); #clean up all old files
	$self->dispatch('install');
    }
}

#replacement for Module::Build's ACTION_install
sub ACTION_install {
    my $self = shift;

    unless($self->config_data('build_done')){
	$self->depends_on('build');
    }

    $self->log_info("Installing " . $self->dist_name . "...\n");
    $self->SUPER::ACTION_install();

    my $blib = $self->blib();
    my $pdir = $self->base_dir."/../perl";
    my @files = grep {-f $_} <$blib/config-*>;
    foreach my $file (@files){
	my ($name) = $file =~ /([^\/]+)$/;
	ExtUtils::Install::pm_to_blib({$file => "$pdir/$name"}, $pdir);
    }

    $self->config_data('install_done' => 1);
}

#replaces Module::Build's ACTION_installdeps
#so prereqs install locally inside of .../perl/lib
sub ACTION_installdeps{
    my $self = shift;

    my $prereq = $self->prereq_failures();
    if($prereq && ($prereq->{requires} || $prereq->{recommends})){
	my @perl = map {keys %{$_->{requires}}} $self->prereq_failures();
	my $access = 1 if(-w $Config{installsitelib} && -w $Config{installarchlib});
	if(!$access){
	    my ($root_id, $grp_id) = (getpwnam('root'))[2,3];
	    $access = 1 if($< == $root_id || $> == $root_id);
	}

	my $local;
	if(! $access){
	    $local = $self->y_n("You do not have write access to install missing Modules.\n".
				"I can try and install these locally (i.e. only for package)\n".
				"in the .../<package>/perl/lib directory, or you can run\n".
				"'./Build installdeps' as root or using sudo and try again.\n".
				"Do want to try and build a local installation?", 'N');

	}	

	if(!$access && !$local){
	    print "\n\nWARNING: You do not appear to have write access to install missing\n".
		  "Modules. Please run './Build installdeps' as root or using sudo.\n\n";

	    my $go = $self->y_n("Do you want to continue anyway?", 'N');

	    exit(0) unless($go);
	}

	my %errors;
	foreach my $m (keys %{$prereq->{build_requires}}){    
	    $self->cpan_install($m, $local);
	}
	foreach my $m (keys %{$prereq->{requires}}){    
	    $self->cpan_install($m, $local);
	}

	print "Checking optional dependencies:\n";
	foreach my $m (keys %{$prereq->{recommends}}){    
	    my $go = $self->y_n("Install $m?", 'Y');
	    if($go){
		$errors{$m}++ unless($self->cpan_install($m, $local));
	    }
	}
	
	print "\nRechecking dependencies to see if installation was successful\n";
	$self->check_prereq;
	
	if($self->prereq_failures() && keys %errors){
	    if(! keys %{$prereq->{requires}}){
		print "\nWARNING: You have all required dependencies, but some optional components\n".
		    "failed to install.\n";
	    }
	    else{
		my ($usr_id) = (getpwnam('root'))[2];
		print "\nWARNING: Installation failed (please review any previous errors).\n";
		print "Try installing the missing packages as 'root' or using sudo.\n" if($< != $usr_id);
		print "You may need to configure and install these packages manually.\n";
		return 0;
	    }
	}
    }

    $self->status;
}

#these install individual external algorithms
sub ACTION_samtools{ shift->_exe_action('samtools'); }
sub ACTION_tabix{ shift->_exe_action('tabix'); }
sub ACTION_vcftools{ shift->_exe_action('vcftools'); }
sub ACTION_hg19{ shift->_exe_action('hg19', 'hg19_random.fa'); }
sub ACTION_mpich2{ shift->_exe_action('mpich2', 'mpich2version'); }

#runs all the algorithm installs that are missing
sub ACTION_installexes{
    my $self = shift;

    my $exe_failures = $self->exe_failures();

    return if(! $exe_failures || ! $exe_failures->{exe_requires});
    foreach my $name (keys %{$exe_failures->{exe_requires}}){
	next if(!$exe_failures->{exe_requires}{$name});
	$name = lc($name); #all dispatch parameters are lower case
	$self->dispatch($name);
	$exe_failures = $self->exe_failures(); #solve recusive install
    }
    
    $self->check_exes;
}

#prints out a simle configuration status message
sub ACTION_status {
    my $self = shift;
    
    $self->status;
}

sub ACTION_clean {
    my $self = shift;
    
    foreach my $f (keys %{$self->PL_files}){
	$f =~ s/Build.PL$/Build/;
	next unless(-f $f);
	print "$f...\n";
	$self->do_system("$f clean");
    }
    print "main package...\n";
    $self->SUPER::ACTION_clean();
}

sub ACTION_realclean {
    my $self = shift;

    foreach my $f (keys %{$self->PL_files}){
	$f =~ s/Build.PL$/Build/;
	next unless(-f $f);
	print "$f...\n";
	$self->do_system("$f realclean");
    }
    $self->SUPER::ACTION_realclean();
}

#checks external algorithms to see if they're present. Anologous to check_prereqs
sub check_exes{
    my $self = shift;

    my $exe_failures = $self->exe_failures();
    return if(! $exe_failures);

    print "Checking external program dependencies...\n";
    while(my $cat = each %$exe_failures){
	my $s = ($cat eq 'exe_requires') ? '!' : '*';
	$cat =~ /exe_(.*)/;

	print "  $1:\n";

	while(my $name = each %{$exe_failures->{$cat}}){
	    print "    $s  ".$exe_failures->{$cat}{$name}{message} ."\n";
	}
    }
    print "\nERRORS/WARNINGS FOUND IN PREREQUISITES.  You may wish to install the programs\n".
	"indicated above before proceeding with this installation.\n".
	"Run 'Build installexes' to install missing prerequisites.\n\n";
}

#checks external libraries to see if they're present. Anologous to check_prereqs
sub check_libs{
    my $self = shift;

    my $lib_failures = $self->lib_failures();
    return if(! $lib_failures);

    print "Checking external library dependencies...\n";
    while(my $cat = each %$lib_failures){
	my $s = ($cat eq 'lib_requires') ? '!' : '*';
	$cat =~ /lib_(.*)/;

	print "  $1:\n";

	while(my $name = each %{$lib_failures->{$cat}}){
	    print "    $s  ".$lib_failures->{$cat}{$name}{message} ."\n";
	}
    }
    print "\nERRORS/WARNINGS FOUND IN PREREQUISITES.  You may wish to install the libraries\n".
	"indicated above before proceeding with this installation\n".
	"\nInstaller cannot do this for you.  You will have to do it manually.\n";
}

#returns missing exes, anologous to prereq_failures
{
my $set_PATH = 0;
sub exe_failures {
    my $self = shift;
    my $other = shift || {};
    my %exes = (%{$self->exe_requires}, %{$other});

    if(! $set_PATH && $self->config_data('PATH')){
	$ENV{PATH} = ($ENV{PATH}) ?
	    $ENV{PATH}.":".$self->config_data('PATH') : $self->config_data('PATH');
	$set_PATH = 1;
    }

    my %exe_failures;
    while (my $name = each %exes){
	my @cmds = map {$_ =~ s/\s+$|^\s+//g; $_} split(/,/, $exes{$name});
	my $base = $self->install_destination('exe');
	my $dest = (-d "$base/$name") ? "$base/$name" : "$base/".lc($name);

	my $loc;
	foreach my $cmd (@cmds){
	    last if(($loc) = grep {$_ && -f $_} $self->config_data($cmd));
	    last if(($loc) = grep {-f $_} ("$dest/$cmd", "$dest/bin/$cmd", File::Which::which($cmd)));
	    last if(($loc) = grep {-f $_ && -x $_} ($cmd));
	}

        if(! $loc){
	    $exe_failures{'exe_requires'}{$name}{have} = '<none>';
	    $exe_failures{'exe_requires'}{$name}{message} = "$name is not installed";
	    $exe_failures{'exe_requires'}{$name}{need} = $exes{$name};
	}
    }

    return (keys %exe_failures) ? \%exe_failures : undef;
}
}

#returns missing libs, anologous to prereq_failures
sub lib_failures {
    my $self = shift;
    my $other = shift || {};
    my %libs = (%{$self->lib_requires}, %{$other});

    my $user_extra = $self->config_data('include_extra');
    push(@DynaLoader::dl_library_path, @{$user_extra}) if($user_extra);
    push(@DynaLoader::dl_library_path, @{$self->include_dirs});
    push(@DynaLoader::dl_library_path, '/usr/X11/include', '/usr/X11/include', '/sw/lib', '/sw/include'); #for libpng

    my %lib_failures;
    while (my $name = each %libs){
	my @cmds = map {$_ =~ s/\s+$|^\s+//g; $_} split(/,/, $libs{$name});

	my $loc;
	foreach my $cmd (@cmds){
	    last if(($loc) = grep {-f $_} (DynaLoader::dl_findfile($cmd), $cmd));
	}

        if(! $loc){
	    $lib_failures{'lib_requires'}{$name}{have} = '<none>';
	    $lib_failures{'lib_requires'}{$name}{message} = "$name is not installed";
	    $lib_failures{'lib_requires'}{$name}{need} = $libs{$name};
	}
    }

    return (keys %lib_failures) ? \%lib_failures : undef;
}


#hidden connection entry between ACTION_??? algorithm install
#checks to see if already installed and then run the install method
sub _exe_action{
    my $self = shift;
    my $label = shift;
    my $script = shift;

    if(! $script && ! $self->exe_requires->{$label}){
	$script = $label;
    }

    my $fail = ($script) ? $self->exe_failures({$label => $script}) : $self->exe_failures();
    my @list = keys %{$fail->{exe_requires}} if($fail->{exe_requires});
    if(grep {/^$label$/i} @list){
	$self->_install_exe($label);
    }
    else{
	my $go = $self->y_n("WARNING: $label was already found on this system.\n".
			    "Do you still want me to install $label for you?", 'N');
	$self->_install_exe($label) if($go);
    }
}

#does actual installation of all external algorithms
sub _install_exe {
    my $self = shift;
    my $exe  = shift;
    my $base = $self->install_destination('exe');
    my $path = "$base/$exe";

    #get OS and architecture
    my %os_ok = (Linux_x86_64  => 1,
		 Linux_i386    => 1,
		 Darwin_i386   => 1,
		 Darwin_x86_64 => 1,
		 src           => 1); 
    my ($OS, $ARC) = (POSIX::uname())[0,4];

    #set all pentium achitectures to use i386 (safest choice)
    if($ARC =~ /^i.86$/){
	$ARC = 'i386';
    }

    ($OS, $ARC) = ('src', '') if(! $os_ok{"$OS\_$ARC"}); #use source code for unknown architectures

    #add fink paths just in case
    if($OS eq 'Darwin'){
	$ENV{C_INCLUDE_PATH} .= ":" if($ENV{C_INCLUDE_PATH});
	$ENV{LIBRARY_PATH} .= ":" if($ENV{LIBRARY_PATH});
	$ENV{C_INCLUDE_PATH} .= "/sw/include:/usr/X11/include";
	$ENV{LIBRARY_PATH} .= "/sw/lib:/usr/X11/lib";
    }

    #get url for exectuable to be installed
    my $data;
    open(LOC, '<', $self->base_dir()."/locations")
	or die "ERROR: Could not open locations to download dependencies\n";
    my $line = <LOC>;
    if($line =~ /^\#\#PACKAGE/){
	$data = join('', <LOC>);
	eval $data;
    }
    close(LOC);

    #prerequisite installation directory
    if(! -d $base){
	mkdir($base) ||
	    die "ERROR could not create directory $base for installing external program dependencies\n";
    }

    #install
    chdir($base);
    my @unlink;
    if($exe eq 'vcftools'){
	#tabix
	my $fail = $self->exe_failures({tabix => 'tabix'});
	$self->_install_exe('tabix') if($fail->{exe_requires}{tabix});

	#vcftools
	File::Path::rmtree($path);
	my $file = "$base/$exe.tar.gz"; #file to save to
        my $url = $data->{$exe}{"$OS\_$ARC"}; #url to blast for OS
	print "Downloading $exe...\n";
        $self->getstore($url, $file) or return $self->fail($exe, $path);
	print "Unpacking $exe tarball...\n";
        $self->extract_archive($file) or return $self->fail($exe, $path);
        push (@unlink, $file);
	my ($dir) = grep {-d $_} <vcftools_*>;
	if(-d $dir){ #this is the source code and must be compiled
	    chdir($dir);
	    print "Configuring $exe...\n";
	    $self->do_system("make") or return $self->fail($exe, $path);
	}
        chdir($base);
	File::Copy::move($dir, $exe) or return $self->fail($exe, $path);
	return $self->fail($exe, $path) if(! -f "$path/bin/vcftools");

	#change back to proper directory or config_data doesn't work
	chdir($self->base_dir());
	$self->config_data(vcftools => "$path/bin/vcftools");
    }
    elsif($exe eq 'tabix'){
	#tabix
	File::Path::rmtree($path);
	my $file = "$base/$exe.tar.bz2"; #file to save to
        my $url = $data->{$exe}{"$OS\_$ARC"}; #url to blast for OS
	print "Downloading $exe...\n";
        $self->getstore($url, $file) or return $self->fail($exe, $path);
	print "Unpacking $exe tarball...\n";
        $self->extract_archive($file) or return $self->fail($exe, $path);
        push (@unlink, $file);
	my ($dir) = grep {-d $_} <tabix-*>;
	if(-d $dir){ #this is the source code and must be compiled
	    chdir($dir);
	    print "Configuring $exe...\n";
	    $self->do_system("make") or return $self->fail($exe, $path);
	}
        chdir($base);
	File::Copy::move($dir, $exe) or return $self->fail($exe, $path);;
	return $self->fail($exe, $path) if(! -f "$path/tabix");

	#change back to proper directory or config_data doesn't work
	chdir($self->base_dir());
	$self->config_data(tabix => "$path/tabix");
    }
    elsif($exe eq 'samtools'){
	#samtools
	File::Path::rmtree($path);
	my $file = "$base/$exe.tar.bz2"; #file to save to
        my $url = $data->{$exe}{"$OS\_$ARC"}; #url to blast for OS
	print "Downloading $exe...\n";
        $self->getstore($url, $file) or return $self->fail($exe, $path);
	print "Unpacking $exe tarball...\n";
        $self->extract_archive($file) or return $self->fail($exe, $path);
        push (@unlink, $file);
	my ($dir) = grep {-d $_} <samtools-*>;
	if(-d $dir){ #this is the source code and must be compiled
	    chdir($dir);
	    print "Configuring $exe...\n";
	    $self->do_system("make dylib") or return $self->fail($exe, $path);
	    $self->do_system("make") or return $self->fail($exe, $path);
	}
        chdir($base);
	File::Copy::move($dir, $exe) or return $self->fail($exe, $path);;
	return $self->fail($exe, $path) if(! -f "$path/samtools");

	#change back to proper directory or config_data doesn't work
	chdir($self->base_dir());
	$self->config_data(samtools => "$path/samtools");
	$self->config_data(sam_include => $path);
	$self->config_data(sam_lib => $path);
	$self->config_data(include_extra => [$path]);
    }
    elsif($exe eq 'mpich2'){
	#MPICH2
	&File::Path::rmtree($path);
	my $file = "$base/$exe.tar.gz"; #file to save to
	my $url = $data->{$exe}{"$OS\_$ARC"}; #url to blast for OS
	print "Downloading $exe...\n";
	$self->getstore($url, $file) or return $self->fail($exe, $path);
	print "Unpacking $exe tarball...\n";
	$self->extract_archive($file) or return $self->fail($exe, $path);
	push (@unlink, $file);
	my ($dir) = grep {-d $_} <mpich2*>;
	print "Configuring $exe...\n";
	chdir($dir);
	my %shared = (Linux  => '--enable-sharedlibs=gcc',
		      Darwin => '--enable-sharedlibs=osx-gcc --disable-f77 --disable-fc',
		      src    => '');
	$self->do_system("./configure --prefix=$path --enable-shared $shared{$OS}") or return $self->fail($exe, $path);
	$self->do_system("make") or return $self->fail($exe, $path);
	$self->do_system("make install") or return $self->fail($exe, $path);
	chdir($base);
	File::Path::rmtree($dir) or return $self->fail($exe, $path);
    }
    elsif($exe eq 'hg19'){
	#use different install location for data
	chdir($self->base_dir());
	my $base = $self->install_destination('data');
	chdir($base);

	#UCSC hg19
	my $url = $data->{$exe}{"$OS\_$ARC"}; #url to blast for OS
	my $file = "$base/$exe\_random.fa"; #file to save to
	print "Downloading $exe...\n";
        $self->get_ucsc_file($url, $file) or return $self->fail($exe, $file);
	return $self->fail($exe, $file) if(! -f $file);

	#change back to proper directory or config_data doesn't work
	chdir($self->base_dir());
	$self->config_data(hg19 => $file);
    }
    else{
	die "ERROR: No install method defined for $exe in INSTALL::Build::_install_exe.\n";
    }

    #report status
    print "Finished installing $exe.\n";

    #remove all tarballs/etc
    map{unlink($_)} @unlink;

    #change back to proper directory
    chdir($self->base_dir());
}

#for downloading UCSC datafiles to a location
sub get_ucsc_file {
    my $self = shift;
    my $url  = shift;
    my $file = shift; #final destination

    #get names for destination files
    (my $dir = $file) =~ s/[^\/]+$//;
    my ($name) = $url =~ /([^\/]+)$/;
    my $zfile = "$dir/$name";
    my $tag_file = "$dir/.cache";

    #get headers from download URL
    my $ua = LWP::UserAgent->new();
    my $r = $ua->head($url);
    if($r->status_line =~ /^500/){
        warn "WARNING: Failed to connect to UCSC. Please check your network.\n";
	return 0;
    }

    #process etag header with cache (see if file is different)
    my %cache;
    eval "require AnyDBM_File";
    eval {	
        mkdir($dir) unless(-d $dir);
        tie(%cache, 'AnyDBM_File', $tag_file, O_CREAT|O_RDWR, 0644)
            or die "Can't open cache file $tag_file: $!";
    };
    if($@){ #second try on failure
        unlink $tag_file;
        tie(%cache, 'AnyDBM_File', $tag_file, O_CREAT|O_RDWR, 0644)
            or die "Can't open cache file $tag_file: $!";
    }
    my $etag   = $r->header('etag') || ''; #md5 checksum
    my $ebytes = $r->header('content-length') || ''; #expected file size
    unlink($zfile) if(!defined($cache{"$url|etag"}) || $cache{"$url|etag"} ne $etag);
    unlink($zfile) if(!defined($cache{"$url|ebytes"}) || $cache{"$url|ebytes"} ne $ebytes);
    $cache{"$url|etag"} = $etag;
    $cache{"$url|ebytes"} = $ebytes;
    untie %cache;

    #check file size for existing downloads
    my $bytes = -s $zfile || 0;
    if($ebytes == $bytes){
        return 1 if(-f $file);
    }

    #prepare handlers for data chunks
    my ($FH, $CH, @param);
    if($url =~ /\.gz$/){ #uncrompresses simultanously with download
        @param = (':content_cb' => sub {print $FH $_[0]; print $CH $_[0]});
        if($url =~ /\.tar\.gz$/){
            open($CH, "| gunzip -c - | tar -O -xf - > $file");
        }
        else{
            open($CH, "| gunzip -c - > $file");
        }

        #restart extraction and partial download from given offset
        if($bytes){
            open(my $IN, "< $zfile");
            print $CH $_ while(<$IN>);
	    close($IN);
            push(@param, 'Range' => "bytes=$bytes-");
        }
    }
    else{ #just get file
        @param = (':content_cb' => sub {print $FH shift @_;});
        push(@param, 'Range' => "bytes=$bytes-") if($bytes); #restart partial
    }

    #download the file
    open($FH, ">> $zfile") or die "ERROR: $!";
    $ua->show_progress(1); #show % complete
    $r = $ua->get($url, @param);
    close($FH);
    close($CH) if($url =~ /\.gz$/);

    #check status of download and return location
    if ($r->status_line =~ /^(200|206|416)/) {
        return 1;
    }
    else{
        warn "WARNING: failed on download with status: ".$r->status_line."\n";
	return 0;
    }
}

#fail/cleanup method for installing exes
sub fail {
    my $self = shift;
    my $exe = shift;
    my $path = shift;

    print "\n\nERROR: Failed installing $exe, now cleaning installation path...\n".
	"You may need to install $exe manually.\n\n";

    if(-f $path){
	unlink($path);
    }
    else{
	File::Path::rmtree($path);
    }
}

# install an external module using CPAN prior to testing and installation
# borrowed and modified from BioPerl, has flag for local install
sub cpan_install {
    my ($self, $desired, $local) = @_;

    unless(File::Which::which('make')){
	die "\n\n".
            "ERROR: Cannot find 'make' on your system. If this is a Mac you will need\n".
	    "to install Xcode developer tools before you can continue. This can be\n".
	    "done from the OS X installation disk or downloaded from the App Store.\n".
            "On some systems you must also select 'install command line tools' under\n".
	    "the 'Downloads' section of the 'Xcode >> Preferences' menu\n\n";
    }

    if(! $local){
	my $loc = $self->module_overide($desired);

	if($loc){
	    print "\n\nWARNING: There is another version of the module $desired\n".
		  "installed on this machine that will supercede a globally installed version.\n".
		  "I can only continue by installing a local package-only version. If you want a\n".
		  "global installation of the module, you will have to quit and delete the\n".
		  "offending module.\n".
                  "Location: $loc\n";

	    $local = $self->y_n("Do you want to continue with a local installation?", 'N');

	    die "\nWARNING: You will need to delete $loc before continuing.\n"if(! $local);
	}
    }

    #set up PERL5LIB environmental varable since CPAN doesn't see my 'use lib'
    if($local){
	my $PERL5LIB = $ENV{PERL5LIB} || '';
	$PERL5LIB = $self->base_dir."/../perl/lib:$PERL5LIB";
	$ENV{PERL5LIB} = $PERL5LIB;
    }
    
    #possible temporary install base location for any CPAN requirements
    my $base = $self->base_dir;
    if(-d "$base/inc/perl"){
	unshift(@INC, "$base/inc/perl/lib");
	my $PERL5LIB = $ENV{PERL5LIB} || '';
	$PERL5LIB = "$base/inc/perl/lib:$PERL5LIB";
	$ENV{PERL5LIB} = $PERL5LIB;
    }

    #if CPAN is too old, then install a newer one first
    if(!$self->check_installed_status('CPAN', '1.82')->{ok}){
	if(! -d "$base/inc/perl"){
	    mkdir("$base/inc/perl");
	    unshift(@INC, "$base/inc/perl/lib");
	    my $PERL5LIB = $ENV{PERL5LIB} || '';
	    $PERL5LIB = "$base/inc/perl/lib:$PERL5LIB";
	    $ENV{PERL5LIB} = $PERL5LIB;
	}

	#fill out environment variables
	my $perl_mb = $ENV{PERL_MB_OPT}; #backup
	$ENV{PERL_MB_OPT} = "--install_base $base/inc/perl/ --installdirs site".
            " --install_path libdoc=$base/inc/perl/man --install_path bindoc=$base/inc/perl/man".
            " --install_path lib=$base/inc/perl/lib --install_path arch=$base/inc/perl/lib".
            " --install_path bin=$base/inc/perl/lib/bin --install_path script=$base/inc/perl/lib/bin ";;
	my $perl_mm = $ENV{PERL_MM_OPT}; #backup
	$ENV{PERL_MM_OPT} = "DESTDIR=$base/inc/perl/ INSTALLDIRS=site INSTALLSITEMAN1DIR=man INSTALLSITEMAN3DIR=man".
            " INSTALLSITEARCH=lib INSTALLSITELIB=lib INSTALLSITEBIN=lib/bin INSTALLSITESCRIPT=lib/bin";
	my $prefer = $ENV{PERL_AUTOINSTALL_PREFER_CPAN}; #backup
	$ENV{PERL_AUTOINSTALL_PREFER_CPAN} = 1;
	my $mm_def = $ENV{PERL_MM_USE_DEFAULT}; #backup
	$ENV{PERL_MM_USE_DEFAULT} = 1;

	#CPAN config from local::lib's Makefile.PL
	my $cpan_config_command =
            'my $done; require ExtUtils::MakeMaker;
             my $orig = ExtUtils::MakeMaker->can("prompt");
             *ExtUtils::MakeMaker::prompt = sub ($;$) {
               if (!$done && $_[0] =~ /manual configuration/) {
                 $done++;
                 return "no";
               }
               return $orig->(@_);
             };
             $CPAN::Config->{prefer_installer} = "EUMM";
             CPAN::Config->load;
             unless ($done || -w $CPAN::Config->{keep_source_where}) {
               my $save = $CPAN::Config->{urllist};
               delete @{$CPAN::Config}{keys %$CPAN::Config};
               $CPAN::Config->{urllist} = $save;
               CPAN::Config->init;
             }';
	my $cpan_command = '';
	$cpan_command .= 'force("notest","install","ExtUtils::MakeMaker"); '
	    if(!$self->check_installed_status('ExtUtils::MakeMaker', '6.31')->{ok});
	$cpan_command .= 'force("notest","install","ExtUtils::Install"); '
	    if(!$self->check_installed_status('ExtUtils::Install', '1.43')->{ok});
	$cpan_command .= 'force("notest","install","CPAN"); ';

	#run CPAN via system call
	system($^X, '-MCPAN', '-e', $cpan_config_command);
	system($^X, '-MCPAN', '-e', $cpan_command);

	$ENV{PERL_MB_OPT} = $perl_mb; #restore
	$ENV{PERL_MM_OPT} = $perl_mm; #restore
	$ENV{PERL_MM_USE_DEFAULT} = $mm_def; #restore
	$ENV{PERL_AUTOINSTALL_PREFER_CPAN} = $prefer; #restore
    }

    # Here we use CPAN to actually install the desired module
    require CPAN;
    import MyModule;

    # Save this because CPAN will chdir all over the place.
    my $cwd = getcwd();

    #set up a non-global local module library for package
    my %bak;
    if($local){
	CPAN::HandleConfig->load;
	%bak = (makepl_arg => $CPAN::Config->{makepl_arg},
		mbuildpl_arg => $CPAN::Config->{mbuildpl_arg});
	$CPAN::Config->{makepl_arg} = "DESTDIR=$base/../perl/ INSTALLDIRS=site INSTALLSITEMAN1DIR=man INSTALLSITEMAN3DIR=man".
	    " INSTALLSITEARCH=lib INSTALLSITELIB=lib INSTALLSITEBIN=lib/bin INSTALLSITESCRIPT=lib/bin";
	$CPAN::Config->{mbuildpl_arg} = "--install_base $base/../perl/ --installdirs site".
	    " --install_path libdoc=$base/../perl/man --install_path bindoc=$base/../perl/man".
	    " --install_path lib=$base/../perl/lib --install_path  arch=$base/../perl/lib".
	    " --install_path bin=$base/../perl/lib/bin --install_path script=$base/../perl/lib/bin ";
	$CPAN::Config->{prefs_dir} = "$ENV{HOME}/.cpan/prefs" if(! -w $CPAN::Config->{prefs_dir});
	CPAN::Shell::setup_output();
	CPAN::Index->reload;
    }
    else{
	CPAN::HandleConfig->load;
	%bak = (makepl_arg => $CPAN::Config->{makepl_arg},
		mbuildpl_arg => $CPAN::Config->{mbuildpl_arg});
	$CPAN::Config->{makepl_arg} = "INSTALLDIRS=site";
	$CPAN::Config->{mbuildpl_arg} = "--installdirs site";
	CPAN::Shell::setup_output();
	CPAN::Index->reload;
    }

    #install YAML if needed to avoid other installation issues with prereqs
    CPAN::Shell->force('notest', 'install', 'YAML') if (! $self->check_installed_status('YAML', '0')->{ok});

    #CPAN::Shell->expand("Module", $desired)->cpan_version <= 2.16;
    #CPAN::Shell->install($desired);
    CPAN::Shell->force('notest', 'install', $desired);

    #restore old CPAN settings
    $CPAN::Config->{makepl_arg} = $bak{makepl_arg};
    $CPAN::Config->{mbuildpl_arg} = $bak{mbuildpl_arg};
    CPAN::Shell::setup_output();

    my $ok;
    my $expanded = CPAN::Shell->expand("Module", $desired);
    if ($expanded && $expanded->uptodate) {
	print "$desired installed successfully\n";
	$ok = 1;
    }
    else {
	print "$desired failed to install\n";
	$ok = 0;
    }
    
    chdir $cwd or die "Cannot chdir() back to $cwd: $!";
    return $ok;
}

#untars a package. Tries to use tar first then moves to the perl package untar Module.
sub extract_archive {
    my $self = shift;
    my $file = shift;

    return 0 if(! $file);
    
    if(File::Which::which('tar')){
	my $command;
	my $u = scalar getpwuid($>);
	my $g = scalar getgrgid($));
	if($file =~ /\.gz$|\.tgz$/){
	    $command = "tar -zxm -f $file";
	}
	elsif($file =~ /\.bz2?$|\.tbz2?$/){
	    $command = "tar -jxm -f $file";
	}
	else{
	    $command = "tar -xm -f $file";
	}
	$command .= " --owner $u --group $g" unless((POSIX::uname())[0] =~ /darwin/i);

	return $self->do_system($command); #fast
    }
    else{
	die "ERROR: Archive::Tar required to unpack missing executables.\n".
	    "Try running ./Build installdeps first.\n\n"
	    if(!$self->check_installed_status('Archive::Tar', '0')->{ok});

	return (Archive::Tar->extract_archive($file)) ? 1 : 0; #slow
    }
}

#downloads files from the internet.  Tries to use wget, then curl,
#and finally LWP::Simple
sub getstore {
    my $self = shift;
    my $url = shift;
    my $file = shift;
    my $user = shift;
    my $pass = shift;

    if(File::Which::which('wget')){ #Linux
	my $command = "wget $url -c -O $file --no-check-certificate";
	$command .= " --user $user --password $pass" if(defined($user) && defined($pass));
	return $self->do_system($command); #gives status and can continue partial
    }
    elsif(File::Which::which('curl')){ #Mac
	my $command = "curl --connect-timeout 30 -f -L $url -o $file";
	$command .= " --user $user:$pass" if(defined($user) && defined($pass));
	my $continue = " -C -";

	#gives status and can continue partial
	my $stat = $self->do_system($command . $continue);
	#just redo if continue fails
	$stat = $self->do_system($command) if(! $stat);
	return $stat;
    }
    else{
	die "ERROR: LWP::Simple required to download missing executables\n".
	    "Try running ./Build installdeps first.\n\n"
	    if(!$self->check_installed_status('LWP::Simple', '0')->{ok});

	$url =~ s/^([^\:]\;\/\/)/$1\:\/\/$user\:$pass\@/ if(defined($user) && defined($pass));
	return LWP::Simple::getstore($url, $file); #just gets the file with no features
    }
}

#prints a nice status message for package configuration and install
sub status {
    my $self = shift;

    my @perl = map {keys %{$_->{requires}}} $self->prereq_failures();
    my @exes = map {keys %{$_->{exe_requires}}} $self->exe_failures();
    my @libs = map {keys %{$_->{lib_requires}}} $self->lib_failures();

    my $dist_name = $self->dist_name;
    my $dist_version = $self->dist_version;

    my $mpi = ($self->feature('mpi_support')) ? 'ENABLED' : 'DISABLED';
    my $stat = 'CONFIGURATION OK';
    $stat = 'MISSING PREREQUISITES' if(@perl || @exes || @libs);
    $stat = 'INSTALLED' if($self->config_data('install_done'));

    print "\n\n";
    print "==============================================================================\n";
    print "STATUS $dist_name $dist_version\n";
    print "==============================================================================\n";
    print "PERL Dependencies:\t";
    print ((@perl) ? 'MISSING' : 'VERIFIED');
    print"\n";
    print "\t\t  !  ". join("\n\t\t  !  ", @perl) ."\n\n" if(@perl);
    print "External Programs:\t";
    print ((@exes) ? 'MISSING' : 'VERIFIED');
    print "\n";
    print "\t\t  !  ". join("\n\t\t  !  ", @exes) ."\n\n" if(@exes);
    print "External C Libraries:\t";
    print ((@libs) ? 'MISSING' : 'VERIFIED');
    print "\n";
    print "\t\t  !  ". join("\n\t\t  !  ", @libs) ."\n\n" if(@libs);
    #print "MPI SUPPORT:\t\t";
    #print $mpi;
    #print "\n";
    print $self->dist_name." PACKAGE:\t";
    print $stat;
    print "\n";

    print "\n\nImportant Commands:\n".
        "\t./Build installdeps\t\#installs missing PERL dependencies\n".
        "\t./Build installexes\t\#installs all missing external programs\n".
        "\t./Build install\t\t\#installs this package, (".$self->dist_name.")\n".
        "\t./Build status\t\t\#Shows this status menu\n\n".
        "Other Commands:\n".
        "\t./Build samtools\t\#installs Samtools\n".
        "\t./Build vcftools\t\#installs VCF-tools (also installs Tabix)\n".
        "\t./Build tabix\t\t\#installs Tabix\n".
        "\t./Build hg19\t\t\#gets human reference genome from UCSC\n";
        #"\t./Build mpich2\t\t\#installs MPICH2 (but manual install recommended)\n";
}

#test if there is another version of the module overriding the CPAN install
sub module_overide {
    my $self = shift;
    my $desired = shift;

    my $mod = $desired; #holds expected .pm file name
    $mod =~ s/\:\:/\//g;
    $mod .= '.pm';
    
    my $test=  qq(\@INC = qw($Config{installsitelib}
			     $Config{installsitearch}
			     $Config{installvendorlib}
			     $Config{installvendorarch}
			     $Config{installprivlib}
			     $Config{installarchlib});
		  require $desired;
		  print \$INC{\"$mod\"};
		  );
    
    my $ok = `$^X -e 'eval q{$test} or exit 1'`;
    my $loc = $self->module_loc($desired) if($ok);
    
    return ($loc && $loc ne $ok) ? $loc : undef;
}

#gets the location of a module
sub module_loc {
    my $self = shift;
    my $desired = shift;

    return if(! $desired);

    eval "require $desired"; #loads module into \%INC    

    $desired =~ s/\:\:/\//g;
    $desired .= ".pm";

    return $INC{$desired};
}

sub svn_w_args {
    my $self = shift;
    my $param = shift;
    my $o_args = shift;

    my $svn = File::Which::which("svn");
    if($svn){
	#get message off command line
	$svn .= " $param";
	if(defined($o_args)){
	    $svn .= " $o_args";
	}
	else{
	    foreach my $arg (@{$self->args->{ARGV}}){
		if($arg =~ /[\s\t]/){
		    $arg =~ s/\'/\\\'/g;
		    $arg = "'$arg'" 
		    }
		$svn .= " $arg";
	    }
	}
	$svn .= " ".$self->base_dir."/../";

	$self->do_system($svn);
    }
    else{
	die "ERROR: Cannot find the executable svn (subversion respository tool)\n";
    }
}

sub load_w_o_header {
    my $file = shift;

    my $data;
    open(IN, "< $file");
    while(my $line = <IN>){	
	#strip of perl shebang portiion 
	while($line =~ /^\#!.*perl/ ||
	      $line =~ /exec\s+.*perl/ ||
	      $line =~ /if\s*0\;/ ||
	      $line =~ /^[\s\t\n]+$/
	      ){
	    $line = <IN>;
	}
	$data = join('', $line, <IN>);
    }
    close(IN);

    return $data;
}

sub sync_bins {
    my $self = shift;
    my $cwd = $self->base_dir;

    my $BIN  = "$cwd/../bin";
    my $sbin = "$cwd/bin";
    my $ibin = "$cwd/inc/bin";

    my @ifiles = map {$_ =~ /([^\/]+)$/} grep {-f $_ && !/\~$/ && !/\.PL$/} <$ibin/*>;
    my @sfiles = map {$_ =~ /([^\/]+)$/} grep {-f $_ && !/\~$/ && !/\.PL$/} <$sbin/*>;

    foreach my $file (@ifiles, @sfiles){
	my $bfile = "$BIN/$file";
	my $rfile = (-e "$ibin/$file") ? "$ibin/$file" : "$sbin/$file";

	next if(! -e $bfile);
	#-w permission must be set before these files are edited by the user
	next unless(sprintf("%04o", (stat($bfile))[2] & 07777) =~ /[2367]/);

	my $bdata = load_w_o_header($bfile);
	my $rdata = load_w_o_header($rfile);

	#scripts have been altered by user
	if($bdata ne $rdata){
	    my $bmod = (stat($bfile))[9];
	    my $rmod = (stat($rfile))[9];
	    
	    if($bmod > $rmod){
		print "copying $bfile  -->  $rfile\n";

		#backup incase of failure
		File::Copy::move($rfile, "$rfile.bk~");

		if(open(IN, "> $rfile")){
		    print IN "#!/usr/bin/perl\n\n";
		    print IN $bdata;
		    close(IN);
		}

		#restore on failure
		if(! -f $rfile && -f "$rfile.bk~"){
		    File::Copy::move("$rfile.bk~", $rfile);
		}
		else{
		    unlink("$rfile.bk~");
		}
	    }
	}	
    }
}

sub check_update_version {
    my $self = shift;

    #get current subversion version
    my ($svn) = `svn info` =~ /Revision\:\s*(\d+)/;
    die "ERROR: Could not query subversion repository\n" if(!$svn);

    #get old version information for last stable release
    open(IN, "< version") or die "ERROR: Could not open version file\n";
    my $data = join("\n", <IN>);
    my ($old_svn) = $data =~ /\$SVN=(\d+)/;
    my ($old_version) = $data =~ /\$VERSION=([\d\.]+)/;
    close(IN);

    #check if update is really needed
    my $version = $old_version;
    if($old_svn == $svn){
	print "Package is already up to date as stable release $version\n";

	return $version;
    }
    else{
	#set new version
	$old_version =~ /(.*)\.(\d+)$/;
	$version = $1;
	my $s = $2; #sub version
	my $n = sprintf ('%02s', $s + 1); #new sub version
	$s = ".$s"; #add decimal
	$n = ".$n"; #add decimal
	
	#if version iteration results in lower value then make sub iterator
	#this means major version numbers can only be changed by the user
	if($n < $s){
	    $n = "$s.01";
	}
	$version .= $n;

	#output what will be next version to file
	#then another commit will be performed to
	#sync subverion with the release
	my $commit_svn = $svn;
	do{
	    $svn = $commit_svn;
	    $svn++;
	    open(OUT, "> version");
	    print OUT "\$VERSION=$version\n";
	    print OUT "\$SVN=$svn\n";
	    close(OUT);

	    #files to fix version for
	    my $cwd = $self->base_dir;
	    my @files = <$cwd/bin/*.pl>;

	    #changing script version here
	    foreach my $file (@files){
		open(IN, "< $file");
		unlink($file);
		open(OUT, "> $file");
		while(my $line = <IN>){
		    $line =~ s/\$VERSION\s*\=\s*\'[\d\.]+\'/\$VERSION = \'$version\'/;
		    print OUT $line;
		}
		close(OUT);
		close(IN);
	    }

	    $self->svn_w_args('commit', "-m \"Stable release version $version\"");
	    $self->svn_w_args('update', '');
	    ($commit_svn) = `svn info` =~ /Revision\:\s*(\d+)/;
	    die "ERROR: Could not query subversion repository\n" if(!$commit_svn);

	    my $svn_server = 'svn://some_server/package';
	    my $copy_args = "$svn_server/trunk $svn_server/tags/Version_$version\_r$svn";
	    my $copy_message = "Adding tags/Version_$version\_r$svn";
	    $self->svn_w_args('copy', "$copy_args -m '$copy_message'");
	}while($svn != $commit_svn);

	print "Package has been updated to stable release $version\n";

	return $version;
    }
}

sub safe_prompt {
    require Term::ReadKey;

    my $self = shift;
    my $m = shift;
    my $d = shift || '';

    print "$m [".('*'x(length($d)))." ]";

    my $key = 0;
    my $r = "";
    #Start reading the keys
    Term::ReadKey::ReadMode(4); #Disable the control keys (raw mode)

    while(ord($key = Term::ReadKey::ReadKey(0)) != 10) { #This will continue until the Enter key is pressed (decimal value of 10)
	if(ord($key) == 127 || ord($key) == 8) { #DEL/Backspace was pressed
	    if(length($r) > 0){
		#1. Remove the last char from the password
		chop($r);
		#2 move the cursor back by one, print a blank character, move the cursor back by one
		print "\b \b";
	    }
	} elsif(ord($key) < 32) {
	    # Do nothing with these control characters
	} else {
	    $r = $r.$key;
	    print "*";
	}
    }
    print "\n"; #because the user pressed enter
    Term::ReadKey::ReadMode(0); #Reset the terminal once we are done

    $r = $d if(length($r) == 0);
    return $r; #Return the response
}

1;
