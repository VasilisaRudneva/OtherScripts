#! /usr/bin/env perl

use warnings;
use strict; 
use Cwd;
use File::Copy;

my $input_file = $ARGV[0]; 
my $num = $ARGV[1];

open IN, "$input_file" or die "cannot open file: $input_file";
my @array=<IN>; 
close IN;


my @result;
my $pos=0;
while (@result<$num) 
	{
        $pos++ while (rand(@$array-$pos)>($num-@result));
        push @result,$array->[$pos++];
    }
foreach my $this_line (@result)
{print "$this_line";}
