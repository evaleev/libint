#!/usr/bin/perl

my $timencycles = 10;

my $nflops;  my $nints;
my $wtime = 0.0;
while ($wtime < 1.0 || $wtime > 10.0) {
  printf "$timencycles\n";
  my $time_output = `./time_eri $timencycles 2>&1`;
  ($nflops) = ($time_output =~ /nflops = (\d*)/g);
  ($nints)  = ($time_output =~ /nints = (\d*)/g);
  ($wtime) = ($time_output =~ /wtime = ([\d\.e-]*)/g);
  if ($wtime < 0.2) {
    $timencycles *= 10;
    next;
  }
  if ($wtime < 1.0) {
    $timencycles *= 2;
    next;
  }
  if ($wtime > 10.0) {
    $timencycles /= 2;
    next;
  }
}

my $floprate = $nflops / $wtime;
my $costperint = $wtime / ($timencycles * $nints);
printf STDOUT "FLOP rate = %10.3f GFLOPs Cost per integral = %15.3f nanosec\n", $floprate/1000000000, $costperint * 1000000000.0;


