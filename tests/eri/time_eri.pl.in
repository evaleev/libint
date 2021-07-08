#!/usr/bin/perl

#
# Runs standard timing ERI runs
#

use Getopt::Long;

my $compdest = "";
my $rundest = "";
&GetOptions("compdest=s" => \$compdest,
            "rundest=s" => \$rundest);

(usage() and die) if ($#ARGV >= 0);

$default_vecmeth = 0;
$default_veclen = 64;
@jobs = ();
add_job(0,0,1,0);
add_job(0,0,2,0);
add_job(1,0,1,0);
add_job(2,0,2,0);
add_job(3,0,3,0);
add_job(4,0,4,0);
add_job(0,0,1,1);
add_job(0,0,2,2);
add_job(0,0,3,3);
add_job(0,0,4,4);
add_job(1,1,1,1);

add_case_vecmeth(1);
add_case_veclen(1);

my ($comphost,$compdir) = ($compdest =~ /(\w):(\w)/g);
my ($runhost,$rundir) = ($rundest =~ /(\w):(\w)/g);


foreach my $job (@jobs) {
  my $la = $job->{LA};
  my $lb = $job->{LB};
  my $lc = $job->{LC};
  my $ld = $job->{LD};
  my $unroll = 100000;
  my $veclen = $job->{VECLEN};
  my $vecmeth = $job->{VECMETH};
  
  printf STDOUT "$la $lb $lc $ld $veclen $vecmeth\n";
  
  if ($rundest eq "") {
    if ($compdest eq "") {
      system("./test_eri.pl --timencycles=100 $la $lb $lc $ld $unroll $veclen $vecmeth");
    }
    else {
      die "compdest specified but rundest is not";
    }
  }
  else {
    if ($compdest eq "") {
      system("./test_eri.pl $la $lb $lc $ld $unroll $veclen $vecmeth");
    }
    else {
      system("./test_eri.pl --scpexport=compdest $la $lb $lc $ld $unroll $veclen $vecmeth");
      system("ssh $comphost \"cd $compdir; make time_eri\"");
    }
    system("ssh $runhost \"cd $rundir; ./run_time_eri.pl\"");
  }
}

exit(0);

sub usage {

  printf STDERR "USAGE: time_eri.pl [options]\n";
  printf STDERR "       Options:\n";
  printf STDERR "         --compdest=host:dir  -- generate source on this host, transfer data to host:dest, and compile time_eri there\n";
  printf STDERR "         --rundest=host:dir   -- run time_eri on host in directory dir\n";

}

sub copy_job {
  
  my ($job) = @_;
  
  my $new_job = {};
  $new_job->{LA} = $job->{LA};
  $new_job->{LB} = $job->{LB};
  $new_job->{LC} = $job->{LC};
  $new_job->{LD} = $job->{LD};
  $new_job->{VECLEN} = $job->{VECLEN};
  $new_job->{VECMETH} = $job->{VECMETH};
  bless($new_job);
  return $new_job;
}

sub add_job {
  my ($la,$lb,$lc,$ld) = @_;
  
  my $job = {};
  $job->{LA} = $la;
  $job->{LB} = $lb;
  $job->{LC} = $lc;
  $job->{LD} = $ld;
  $job->{VECLEN} = $default_veclen;
  $job->{VECMETH} = $default_vecmeth;
  bless($job);
  push @jobs, $job;
}

sub add_case_vecmeth {
  my ($vecmeth) = @_;
  
  my @orig_jobs = @jobs;
  foreach my $job (@orig_jobs) {
    my $new_job = copy_job($job);
    $new_job->{VECMETH} = $vecmeth;
    push @jobs, $new_job;
  }
}

sub add_case_veclen {
  my ($veclen) = @_;
  
  my @orig_jobs = @jobs;
  foreach my $job (@orig_jobs) {
    my $new_job = copy_job($job);
    $new_job->{VECLEN} = $veclen;
    push @jobs, $new_job;
  }
}

