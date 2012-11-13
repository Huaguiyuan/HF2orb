#!/usr/bin/perl -w
use strict;

my ($beginJob,$endJob,$label)=@ARGV;

for (my $i=$beginJob;$i<$endJob+1;$i++) {
	open(FOUT,">job$label$i.sge") or die "Cannot open job$label$i.sge for writing: $!\n";

print FOUT<<EOF;
#\$ -N HF$label$i
#\$ -q medium*
#\$ -cwd
#\$ -l mem=2G
./hf input$label$i.inp
EOF
close(FOUT);

	system("qsub < job$label$i.sge");
}

