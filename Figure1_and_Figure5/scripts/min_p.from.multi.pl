use strict;
use warnings;

my ($files, $out) = @ARGV;

my %min;

open FL, $files;
while (my $path = <FL>) {
    chomp($path);
    
    open AL, "tail -n +2 $path|";
    while (<AL>) {
        chomp;
        my @tmp = split;
        if (!defined $min{$tmp[0]}{$tmp[2]} or $tmp[11] < $min{$tmp[0]}{$tmp[2]}) {
            $min{$tmp[0]}{$tmp[2]} = $tmp[11];
        }
    };
    close AL;
};
close FL;

open FLS, ">$out";
foreach my $chr (sort {$a<=>$b} keys %min) {
    foreach my $ps (sort {$a<=>$b} keys %{$min{$chr}}) {
        print FLS "si$chr:$ps\t$chr\t$ps\t$min{$chr}{$ps}\n";
    };
};
close FLS;

