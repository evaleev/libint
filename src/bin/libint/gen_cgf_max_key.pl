#!/usr/bin/perl

my $lmax = 20;
my $ncgf = 0;

printf STDOUT "{ 0";
for(my $l=0; $l<=$lmax; $l++) {
  $ncgf += ($l+1)*($l+2)/2;
  printf STDOUT ", $ncgf";
}
printf STDOUT "}\n";

